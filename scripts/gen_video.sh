#!/bin/sh
ffmpeg -r 30 -i ./output/pngs_vis0/%*.png -pix_fmt yuv420p -vcodec libx264 -preset slower ./output/video/vis0.mp4