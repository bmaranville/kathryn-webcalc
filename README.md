**webcalc**

to compile, install the latest version of emsdk_portable (available at http://webassembly.org/getting-started/developers-guide/ )

then, at the command line execute e.g.
emcc --bind -s DISABLE_EXCEPTION_CATCHING=0 SimpleHexagonal7_bbm.cpp -o hex.js -O3
