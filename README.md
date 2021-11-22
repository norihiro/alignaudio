# Audio alignment tool
This is a tool to align audio files.

## Usage example
The primary application is to create a video file from separately recorded video and audio devices.
The steps are as below.
1. Prepare two input files.
   - video file `video-cam.mp4`
   - audio file `audio.wav`
2. Extract audio from the video file.
   ```
   ffmpeg -i video-cam.mp4 video-cam-audio.wav
   ```
3. Align audio and generate aligned audio file.
   ```
   alignaudio video-cam-audio.wav audio.wav -o audio-aligned.wav
   ```
4. Mux video and audio.
   ```
   ffmpeg -i video-cam.mp4 -i audio-aligned.wav -map 0:v -map 1:a -c:v copy video.mp4
   ```
5. Cleanup files.
   ```
   rm video-cam-audio.wav audio-aligned.wav
   ```

## Options
| Option | Description |
| ------ | ----------- |
| `-d align-data-file` | Output results of internal calculation. |
| `-c` | Compensate clock drift |

## Details of the algorithm
This program sweeps audio offset and calculates correlation of the two audio files for each offset.
The sweep will be performed 3 times, coarse sweep at first time, finer sweep, and fine sweep.

The sweep results can be output as a data file readable by GNUplot.
1. Use `-d` option to output the data file.
   ```
   alignaudio -d align.dat ...
   ```
2. Run GNUplot to plot the results of each sweep.
   ```
   set terminal svg
   set output 'align.svg'
   set multiplot layout 3,1
   set xlabel 'time [s]'
   plot 'align.dat' index 0 u 1:2 t 'First sweep'
   set xlabel 'time [s]'
   plot 'align.dat' index 1 u 1:2 t 'Second sweep'
   set xlabel 'time [s]'
   plot 'align.dat' index 2 u 1:2 t 'Third sweep'
   ```

## Clock drift compensation
Since consumer products don't have clock synchronization mechanism,
some tools might suffer drift of the clock source.
This program provides an option to compensate for the clock drift.


## License
This code is released under the General Public License version 3.
See the file 'COPYING' for the full license.
