#!/bin/sh
ffmpeg -r 30 -i ./output/pngs_norigidbody_vis0.02/%*.png -pix_fmt yuv420p -vcodec libx264 -preset slower ./output/video/norigidbody_vis0.02.mp4