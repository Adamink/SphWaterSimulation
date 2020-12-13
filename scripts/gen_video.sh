#!/bin/sh
ffmpeg -r 30 -i ./output/5_frame_pngs/%*.png -pix_fmt yuv420p -vcodec libx264 -preset slower ./output/video/5_frame.mp4