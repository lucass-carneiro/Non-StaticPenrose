# Debug Build instructions

```
mkdir Debug && cd Debug
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-11-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build . -j20
```

# Release Build instructions

```
mkdir Release && cd Release
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-11-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . -j20
```