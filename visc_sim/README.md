Implementation of particle based viscoelastic fluid simulation

Build instructions:

```
mkdir build
cd build
cmake ..
make
visc_sim.exe <scene number> <max num iterations>
```

The code only works on linux.

Executable will generate a directory with png files that can be encoded using ffmpeg or viewed any other way.
