#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

// #define debug(fmt, ...)
#define debug(fmt, ...) fprintf(stderr, fmt, ##__VA_ARGS__)

struct alignaudio_config
{
	const char *file1;
	const char *file2;
	const char *file_out;
	FILE *fp_data;

	int n_average_for_amplitude;
	int n1_average_for_amplitude;
	int n2_average_for_amplitude;

	int n_samples_base_for_amplitude;

	bool compensate_drift;
};

static void alignaudio_config_init(struct alignaudio_config *config)
{
	memset(config, 0, sizeof(struct alignaudio_config));
	config->n_average_for_amplitude = 65536;
	config->n1_average_for_amplitude = 4096;
	config->n2_average_for_amplitude = 256;
	config->n_samples_base_for_amplitude = 48000 * 2 * 60; // 1 minutes
}

static int parse_arguments(struct alignaudio_config *config, int argc, char **argv)
{
	for (int i=1; i<argc; i++) {
		char *ai = argv[i];
		if (*ai=='-') {
			int c;
			while ((c=*++ai)) switch (c) {
				case 'c':
					config->compensate_drift = 1;
					break;
				case 'o':
					if (i+1 >= argc) {
						fprintf(stderr, "Error: -o requires output file name.\n");
						return 29;
					}
					config->file_out = argv[++i];
					break;
				case 'd':
					if (i+1 >= argc) {
						fprintf(stderr, "Error: -o requires output file name.\n");
						return 29;
					}
					config->fp_data = fopen(argv[++i], "w");
					break;
				default:
					fprintf(stderr, "Error: unkown option -%c\n", c);
					return 35;
			}
		}
		else if(!config->file1) {
			config->file1 = ai;
		}
		else if(!config->file2) {
			config->file2 = ai;
		}
	}
	return 0;
}

struct wave_data
{
	uint8_t *data;
	size_t file_length;
	int16_t *raw; // point to the start address in data
	size_t raw_length; // number of samples.

	uint16_t n_channel;
	uint32_t sample_rate;
};

struct alignaudio
{
	struct wave_data data1;
	struct wave_data data2;

	int cand_min, cand_max;
	int center;
	int center_drift1, center_drift2;
};

static inline int wave_comp_sz4(uint8_t *d, char *s) { return memcmp(d, s, 4); }
static inline int wave_read_uint16(uint8_t *d) { return d[0] | (d[1]<<8); }
static inline uint32_t wave_read_uint32(uint8_t *d) { return d[0] | (d[1]<<8) | (d[2]<<16) | (d[3]<<24); }

int load_wave(struct wave_data *w, const char *filename)
{
	FILE *fp = fopen(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t len = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	if (w->data) free(w->data);
	w->data = malloc(len);
	w->file_length = len;
	fread(w->data, 1, len, fp);
	fclose(fp);

	size_t offset = 0;
	bool found_wave = 0;
	bool found_data = 0;
	bool in_riff = 0;
	while (!found_data) {
		uint8_t *p = w->data + offset;
		if (!in_riff && !wave_comp_sz4(p, "RIFF")) {
			int s = wave_read_uint32(w->data + offset + 4);
			debug("load_wave: %s: read RIFF offset=%d size=%d\n", filename, (int)offset, s);
			offset += 8;
			in_riff = 1;
		}
		else if(in_riff && !wave_comp_sz4(p, "WAVE")) {
			debug("load_wave: %s: found WAVE header offset=%d\n", filename, (int)offset);
			offset += 4;
			found_wave = 1;
		}
		else if(in_riff && found_wave && !wave_comp_sz4(p, "fmt ")) {
			int s = wave_read_uint32(p + 4);
			debug("load_wave: %s: found fmt  chunk offset=%d size=%d\n", filename, (int)offset, s);
			offset += 8 + s;
			int format = wave_read_uint16(p + 8); // must be 1: linear PCM
			if (format!=1) {
				fprintf(stderr, "Error: cannot support non linear PCM files. format=%d\n", format);
				return 114;
			}
			w->n_channel = wave_read_uint16(p + 10);
			w->sample_rate = wave_read_uint32(p + 12);
		}
		else if(in_riff && !wave_comp_sz4(p, "LIST")) {
			int s = wave_read_uint32(p + 4);
			debug("load_wave: %s: found LIST chunk offset=%d size=%d\n", filename, (int)offset, s);
			offset += 8 + s;
		}
		else if(found_wave && !wave_comp_sz4(p, "data")) {
			size_t s = wave_read_uint32(p + 4);
			debug("load_wave: %s: found data chunk offset=%d size=%d\n", filename, (int)offset, (int)s);
			w->raw = (int16_t*)(p + 8); // assume the file format matches with type and system endian.
			w->raw_length = s / 2; // assume 16bit
			found_data = 1;
		}
		else {
			fprintf(stderr, "Error: load_wave(%s): unkown chunk %X %X %X %X\n", filename, p[0], p[1], p[2], p[3]);
			return 107;
		}
	}

	return 0;
}

int write_wave(struct wave_data *w, const char *filename)
{
	FILE *fp = fopen(filename, "wb");
	fwrite(w->data, 1, w->file_length, fp);
	fclose(fp);
	return 0;
}

static void mkamplitude(uint16_t *dst, int16_t *src, size_t n_dst, size_t n_amplitude)
{
	// calculate average of the signal for each n_amplitude samples.
	while (n_dst--) {
		uint32_t s = 0;
		for (size_t i=0; i<n_amplitude; i++) {
			int16_t v = *src++;
			if (v>0) s += v;
			else s -= v;
		}
		*dst++ = s / n_amplitude;
	}
}

void calculate_by_amplitude(struct alignaudio *aa, struct alignaudio_config *config)
{
	size_t n_env1 = aa->data1.raw_length / config->n_average_for_amplitude;
	uint16_t *env1 = malloc(sizeof(uint16_t) * n_env1);
	mkamplitude(env1, aa->data1.raw, n_env1, config->n_average_for_amplitude);

	size_t n_env2 = aa->data2.raw_length / config->n_average_for_amplitude;
	uint16_t *env2 = malloc(sizeof(uint16_t) * n_env2);
	mkamplitude(env2, aa->data2.raw, n_env2, config->n_average_for_amplitude);

	debug("calculate_by_amplitude: n_env1=%d n_env2=%d\n", (int)n_env1, (int)n_env2);

	// data1: reference of the time, might be shorter
	// data2: target audio

	struct result_s { int offset; uint64_t sp; } *rr = malloc(sizeof(struct result_s) * (n_env1 + n_env2 + 1));
	int n_rr = 0;

	int offset_start = aa->cand_min / config->n_average_for_amplitude;
	int offset_end   = aa->cand_max / config->n_average_for_amplitude;

	debug("calculate_by_amplitude: comparing... n_env2=%d n_env1=%d offset:[%d:%d]\n", (int)n_env2, (int)n_env1, offset_start, offset_end);

	uint64_t sp_max = 0;
	int offset_ret = 0;
	int64_t d0, d1;
	for (int offset=offset_start; offset<=offset_end; offset++) {
		uint64_t sp = 0;
		for(int i=0; i<n_env1; i++) {
			int j = i - offset;
			if (j < 0) continue;
			if (j >= n_env2) break;
			sp += (uint64_t)env1[i] * env2[j] / 256;
		}

		if (sp > sp_max) {
			d0 = sp - sp_max; // assuming sp_max was the previous data
			offset_ret = offset;
			sp_max = sp;
		}
		else if (offset==offset_ret+1) {
			d1 = sp_max - sp;
		}

		rr[n_rr].offset = offset;
		rr[n_rr].sp = sp;
		n_rr++;
	}

	if (config->n_samples_base_for_amplitude > 0) {
		uint64_t sp_sum = 0;
		int i0 = 0, i1 = 0;
		int m = config->n_samples_base_for_amplitude / config->n_average_for_amplitude;

		sp_max = 0;

		while (i1<n_rr && i1-i0<m)
			sp_sum += rr[i1++].sp;

		for (int i=1; i+1<n_rr; i++) {
			if (i1<n_rr && (i1-i)*2<m)
				sp_sum += rr[i1++].sp;

			if (i1-i0 > m)
				sp_sum -= rr[i0++].sp;

			int64_t sp_avg = sp_sum / (i1-i0);
			int64_t sp = (int64_t)rr[i].sp - sp_avg;

			if (sp > 0 && sp > sp_max) {
				d0 = rr[i].sp - rr[i-1].sp;
				d1 = rr[i].sp - rr[i+1].sp;
				offset_ret = rr[i].offset;
				sp_max = sp;
			}
		}
	}

	// for debug
	if (config->fp_data) {
		for(int i=0; i<n_rr; i++) {
			fprintf(config->fp_data, "%f\t%f\n", (double)rr[i].offset*config->n_average_for_amplitude / (48e3*2), (double)rr[i].sp / (double)sp_max);
		}
		fprintf(config->fp_data, "\n\n");
	}

	// qsort(rr, n_rr, sizeof(struct result_s), (int (*)(const void*, const void*))r_comp);

	int center_adjust;
	center_adjust = config->n_average_for_amplitude * (d0-d1) / (d0+d1) / 2;
	debug("center_adjust=%d d0=%f d1=%f\n", center_adjust, (double)d0, (double)d1);

	aa->center = offset_ret * config->n_average_for_amplitude + (center_adjust&~1);
	aa->cand_min = aa->center - 2 * config->n_average_for_amplitude;
	aa->cand_max = aa->center + 2 * config->n_average_for_amplitude;

	debug("calculate_by_amplitude: cand=[%d:%d] center=%d\n", aa->cand_min, aa->cand_max, aa->center);

	free(env1);
	free(env2);
}

void calculate_drift(struct alignaudio *aa, struct alignaudio_config *config)
{
	const int div = 4 * 2; // 48kHz stereo -> 12kHz mono
	int n1 = aa->data1.raw_length / div;
	int *d1 = malloc(sizeof(int) * n1);
	for (int i=0, j=0, k=0, s=0; i<aa->data1.raw_length; i++) {
		s += aa->data1.raw[i];
		if (++k >= div) {
			d1[j++] = s / div;
			s = 0;
			k = 0;
		}
	}
	int n2 = aa->data2.raw_length / div;
	int *d2 = malloc(sizeof(int) * n2);
	for (int i=0, j=0, k=0, s=0; i<aa->data2.raw_length; i++) {
		s += aa->data2.raw[i];
		if (++k >= div) {
			d2[j++] = s / div;
			s = 0;
			k = 0;
		}
	}

	int offset_begin = aa->cand_min/div;
	int offset_end = aa->cand_max/div;
	int offset_c = (offset_begin+offset_end) / 2;
	double s_xy=0, s_x=0, s_y=0, s_x2=0, s_y2=0, s_w=0;
	const int n_blk = 65536;
	// TODO: maybe this calculation is too fuzy.
	// I should make 1s blocks for example, and find the optimum offset, then use least-square method.
	for (int x=n_blk; x+n_blk<n1; x+=n_blk*32) {
		long long w_max = 0;
		int y = offset_c;
		for (int k=offset_begin-offset_c; k<=offset_end-offset_c; k++) {
			long long w = 0;
			for (int i=x-n_blk/2; i<=x+n_blk/2; i++) {
				int j = i - (k+offset_c);
				if (j<0 || n2<=j) continue;

				w += d1[i] * d2[j];
			}

			if (w<0) w = -w;
			if (w > w_max) {
				w_max = w;
				y = k;
			}
		}
		s_x += x;
		s_y += y;
		s_w += 1;
		s_xy += x * y;
		s_x2 += x * x;
		s_y2 += y * y;
	}

	double a = (s_xy * s_w - s_x * s_y) / (s_w * s_x2 - s_x*s_x);
	double b = (s_x2 * s_y - s_xy * s_x) / (s_w * s_x2 - s_x*s_x);

	printf("calculate_drift detected %.2fppm drift (2nd audio has %s clock)\n", a*1e6, (a<0 ? "slower" : "faster"));
	debug("calculate_drift a=%f b=%f s_w=%f\n", a, b, s_w);
	debug("calculate_drift offset at begin %fs\n", (offset_c+b)*div / 96e3);
	debug("calculate_drift offset at end %fs\n", ((offset_c+b)*div + a*aa->data1.raw_length) / 96e3);

	aa->center = ((offset_c+b)*div);
	aa->center &= ~1;
	aa->center_drift1 = aa->data1.raw_length;
	aa->center_drift2 = aa->data1.raw_length * (1.0+a);
	printf("calculate_drift drift1=%d drift2=%d\n", aa->center_drift1, aa->center_drift2);
	if (aa->center_drift1 > aa->center_drift2) debug("calculate_drift %d samples will be added.\n", aa->center_drift1-aa->center_drift2);
	if (aa->center_drift1 < aa->center_drift2) debug("calculate_drift %d samples will be deleted.\n", aa->center_drift2-aa->center_drift1);

	free(d1);
	free(d2);
}

void overwrite_at_center(struct alignaudio *aa)
{
	int offset = aa->center;
	offset &= ~1; // not to swap stereo channels
	long long drift1 = aa->center_drift1;
	long long drift2 = aa->center_drift2;
	debug("overwrite_at_center offset=%d (%fs) cand was [%fs:%fs]\n", offset, offset/96e3, aa->cand_min/96e3, aa->cand_max/96e3);
	int n_duplicate=0, n_skipped=0;

	for (int i=0, j=-offset; i<aa->data1.raw_length;) {
		if (j>=aa->data2.raw_length) break;
		if (0<=j && j<aa->data2.raw_length)
			aa->data1.raw[i] = aa->data2.raw[j];
		// i*drift2 = j*drift1
		if ((j+offset+1)*drift1 < i*drift2) {
			j++;
			n_skipped++;
		}
		else if ((i+1)*drift2 < (j+offset)*drift1) {
			i++;
			n_duplicate++;
		}
		else {
			i ++;
			j ++;
		}
	}
	debug("overwrite_at_center n_skipped=%d n_duplicate=%d\n", n_skipped, n_duplicate);
}

int main(int argc, char **argv)
{
	int ret;
	struct alignaudio_config config; alignaudio_config_init(&config);
	struct alignaudio aa; memset(&aa, 0, sizeof(aa));

	if ((ret = parse_arguments(&config, argc, argv))) {
		fprintf(stderr, "Error: failed to parse arguments.\n");
		return ret;
	}

	if (!config.file1 || !config.file2) {
		fprintf(stderr, "Error: you have to provide exactly two files.\n");
		return 65;
	}

	ret = load_wave(&aa.data1, config.file1);
	if (ret) return ret;

	ret = load_wave(&aa.data2, config.file2);
	if (ret) return ret;

	aa.cand_min = -(int)aa.data2.raw_length;
	aa.cand_max =  (int)aa.data1.raw_length;
	calculate_by_amplitude(&aa, &config);

	config.n_samples_base_for_amplitude = 0;
	config.n_average_for_amplitude = config.n1_average_for_amplitude;
	calculate_by_amplitude(&aa, &config);

	config.n_average_for_amplitude = config.n2_average_for_amplitude;
	calculate_by_amplitude(&aa, &config);

	// TODO: adjust drift of the clock period

	if (config.compensate_drift) {
		calculate_drift(&aa, &config);
	}

	if (config.file_out) {
		overwrite_at_center(&aa);
		write_wave(&aa.data1, config.file_out);
	}

	return 0;
}
