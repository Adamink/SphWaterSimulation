#!/bin/sh
ffmpeg -r 30 -i ./output/meshed_liquid/c0.5p0.05k10.0s0.6/%*.png -pix_fmt yuv420p -vcodec libx264 -preset slower ./output/video.mp4