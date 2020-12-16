curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install splashsurf

# if errors occurs during installation, like:
#  /usr/bin/ld: cannot find Scrt1.o: No such file or directory
#  /usr/bin/ld: cannot find crti.o: No such file or directory
# use the following commands to fix the problem
# sudo apt-get install libc6-dev

sudo apt install python3-pip
pip3 install -U meshio
