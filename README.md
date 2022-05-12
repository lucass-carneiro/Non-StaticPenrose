# Debug Build instructions

```
mkdir Debug && cd Debug
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-11-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -DGRLENSING_BUILD_TESTING=ON
cmake --build . -j20
```
# Profile Build instructions

```
mkdir Profile && cd Profile
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-11-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Profile 
cmake --build . -j20
```

# Release Build instructions

```
mkdir Release && cd Release
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-11-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . -j20
```

# Executing
```
cd Debug
ln -s ../configs/grlensing_config.yaml ./
```