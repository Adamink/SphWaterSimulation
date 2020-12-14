#!/bin/sh
ffmpeg -r 30 -i ./output/png5_div5/%*.png -pix_fmt yuv420p -vcodec libx264 -preset slower ./output/5_frame.mp4