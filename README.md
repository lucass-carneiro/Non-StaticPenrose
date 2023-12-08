# Instrucoes extras

1. Instalar dependencias basicas
  `python`
  `pip`
  `cmake`
  `mpi` (`openmpi-bin`, `libopenmpi-dev`)

2. Clonar o repositorio

```
git clone https://github.com/lucass-carneiro/Non-StaticPenrose
```

3. Criar o `venv`

```
cd Non-StaticPenrose
python -m venv venv             # Cria o venv
source venv/bin/activate        # Ativa o venv
pip install -r requirements.txt # Instala o conan
deactivate                      # Sai do venv (qdo terminar)
```

4. Compilar o `Release`(#release-build-instructions)

# Debug Build instructions

```
mkdir Debug && cd Debug
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-12-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -DGRLENSING_BUILD_TESTING=ON
cmake --build . -j20
```

# Profile Build instructions

```
mkdir Profile && cd Profile
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-12-release
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Profile 
cmake --build . -j20
```

# Release Build instructions

(dentro do venv)

```
mkdir Release && cd Release
conan install ../conan --remote=conancenter --build missing --profile ../conan/toolchain-gcc-12-release # Necessita do venv
cmake .. -G "Unix Makefiles" -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . -j20
```

# Executing

(dentro de Release)

```
cd Relese
ln -s ../configs/grlensing_config.yaml ./
```
