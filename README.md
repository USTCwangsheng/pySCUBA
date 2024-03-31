pySCUBA
=========

pySCUBA is Python bindings to the SCUBA library. These can be built from the SCUBA source code. 

Getting Started Using pySCUBA
============================

## 1. Requirements

-  Operating System: Linux (Recommended)

### 1.1. Install Dependencies

Before installing, make sure the following dependencies are installed:

-  cmake
-  lapack
-  python3.8

Here is how to install cmake:
```sh
sudo apt install cmake gfortran
```

To install [lapack](https://github.com/Reference-LAPACK/lapack), we suggest build from source code. For example, you can download lapack-3.10.1 and compile it following the instructions in lapack-3.10.1/README.md:

```sh
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lapack ..
   cmake --build . -j --target install
```

We assume that the dynamic libraries liblapacke.so installed in` ~/.local/lapack/` if you follow the previous steps .However , if lapack installed in other path ,your can to change the line in our CMakeLists.txt to your path:
```cmake
# if install lapack in other enviroment ,change the line
set(PythonLibs_LIBRARIES ${YOURPATH}/liblapacke.so)
```
and make sure the file lapacke.h in `/usr/local/include`

Python 3.8 dynamic libraries are also necessary for our project. We assume you have created a Python 3.8 environment using Anaconda. 
```sh
conda create --name pySCUBA python=3.8
conda activate pySCUBA
conda install numpy
``` 
The following files, libpython3.8.so and libpython3.8.so.1.0, should be located in `~/anaconda3/envs/pySCUBA/lib/`. If your python3.8 enviroment is installed in other directory. You have to change the line in our CMakeLists.txt to your path.
```cmake
# if install python3.8 in other enviroment ,change the line
set(PythonLibs_LIBRARIES ${YOURPATH}/libpython3.8.so)
```
### 1.2. Install pySCUBA
Enter our project directory and run the installation program:
```sh
python setup.py develop
```
After a few minutes of compiling, the pySCUBA package is installed

To test if pySCUBA is successfully installed, switch to the `/demo/SCUBASketch/` directory and run
```sh
python SCUBAsketch.py
```
If there are no error messages, it indicates that the pySCUBA has been installed successfully. More demos are in the `/demo/SCUBASketch/` folder.