curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install splashsurf
if any error 
  = note: /usr/bin/ld: cannot find Scrt1.o: No such file or directory
          /usr/bin/ld: cannot find crti.o: No such file or directory
          collect2: error: ld returned 1 exit status
          
sudo apt-get install libc6-dev

sudo apt install meshio-tools /？
pip3 install -U meshio