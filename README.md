# Audio alignment tool
This is a tool to align audio files.

## Usage example
The primary application is to create a video file from separately recorded video and audio devices.
Let's say you have video file `video-cam.mp4` and audio file `audio.wav`,
the steps are as below.
1. Extract audio from the video file.
   `ffmpeg -i video-cam.mp4 video-cam-audio.wav`
2. Align audio and generate aligned audio file.
   `alignaudio video-cam-audio.wav audio.wav -o audio-aligned.wav`
3. Mux video and audio.
   `ffmpeg -i video-cam.mp4 -i audio-aligned.wav -map 0:v -map 1:a -c:v copy video.mp4`
4. Cleanup files.
   `rm video-cam-audio.wav audio-aligned.wav`

## License
This code is released under the General Public License version 3.
See the file 'COPYING' for the full license.
