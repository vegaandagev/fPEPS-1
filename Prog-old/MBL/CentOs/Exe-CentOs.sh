#!/bin/bash


sudo cmake  -DBUILD_DOC=on  -DBUILD_ARPACK_SUPPORT=on   -DBUILD_HDF5_SUPPORT=on  -DBUILD_PYTHON_WRAPPER=on  ~/Downloads/uni10/uni10/

sudo make

sudo make install 


sudo python pyUni10/setup.py install

sudo cp -r /usr/local/uni10/lib/libuni10.* /usr/lib/
sudo cp -r /usr/local/uni10/include/uni10* /usr/include/
 
 
sudo python setup.py install
#sudo python  setup.py install
#Do it directly
#Using chmod u+x scriptname make the script executable.


sudo  CC=/opt/intel/compilers_and_libraries/linux/bin/intel64/icc   CXX=/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc       cmake   -D BUILD_WITH_MKL=True       -DBUILD_PYTHON_WRAPPER=on  -DBUILD_ARPACK_SUPPORT=on   -DBUILD_HDF5_SUPPORT=on   ~/Downloads/uni10-develop/uni10/ 



sudo  CC=/opt/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64/icc    CXX=/opt/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64/icpc      cmake   -D BUILD_WITH_MKL=True       -DBUILD_PYTHON_WRAPPER=on  -DBUILD_ARPACK_SUPPORT=on   -DBUILD_HDF5_SUPPORT=on   ~/Downloads/uni10-develop/uni10/

#add to flags
-I /opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/include 
-L /opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin/ -Wl,-rpath=/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin/







#run
icpc *.cpp  -Wl,-rpath=/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64_lin/   -O2    -std=c++11   -luni10




export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64_lin/:$LD_LIBRAR$
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin/:$LD_LIBRAR$
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl/include:$LD_LIBRAR$


You should add more details about your distribution, for example under Ubuntu the right way to do this is to add a custom .conf file to /etc/ld.so.conf.d, for example

sudo gedit /etc/ld.so.conf.d/randomLibs.conf
inside the file you are supposed to write the complete path to the directory that contains all the libraries that you wish to add to the system, for example

/home/linux/myLocalLibs
remember to add only the path to the dir, not the full path for the file, all the libs inside that path will be automatically indexed.

Save and run sudo ldconfig to update the system with this libs.







sudo   CC=/opt/intel/bin/icc    CXX=/opt/intel/bin/icpc      cmake   -D BUILD_WITH_MKL=True       -DBUILD_PYTHON_WRAPPER=on  -DBUILD_ARPACK_SUPPORT=on   -DBUILD_HDF5_SUPPORT=on   /home/rezahagh/Downloads/uni10-develop/uni10/


-I /opt/intel/mkl/include/

sudo gedit /etc/ld.so.conf.d/randomLibs.conf
/opt/intel/composer_xe_2013.2.146/composer_xe_2013.2.146/compiler/lib/intel64/
/opt/intel/composer_xe_2013.2.146/compiler/lib/intel64/
/opt/intel/mkl/lib/intel64/





I had a similar problem, and resolved it using @Paul's answer as a hint. I needed to use python2.7 to compile an older library, but cmake keeps picking up my python3.2 libraries (and executable).

First, I ran cmake with default options, then edited the CMakeCache.txt file which it generated. I did it this way primarily because I didn't know the proper -D... incantations to cause cmake to get the python library and include paths, etc right in the first place.

In my CmakeCache.txt, I found lines like this

Path to a program
PYTHON_EXECUTABLE:FILEPATH=/usr/bin/python
Path to a directory
PYTHON_INCLUDE_DIR:PATH=/usr/include/python3.2
Path to a library
PYTHON_LIBRARY:FILEPATH=/usr/lib/libpython3.2.so
And replaced every occurrence of python3.2 with python2.7. I also had to rename the PYTHON_EXECUTABLE to use python2.7, since python is a symlink to python3.2 on my system.

Then I reran cmake. Because it prefers its cached values to actually looking for the libraries, this should work in all cases. At least, it did in mine.



Try to add -DPYTHON_EXECUTABLE:FILEPATH=/path/to/python2.7 It might be a path problem?

Also could specify the path to your python library,use your version that you want:

 cmake -DPYTHON_LIBRARIES=/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib .


 cp -r /usr/local/lib/python2.7/site-packages/pyUni* /usr/lib/python2.6/site-packages/
//Path to a program.
PYTHON_EXECUTABLE:FILEPATH=/usr/bin/python2.6

//Path to a file.
PYTHON_INCLUDE_DIR:PATH=/usr/include/python2.6

//Path to a library.
PYTHON_LIBRARY:FILEPATH=/usr/lib/python2.6/



wget --no-check-certificate https://bootstrap.pypa.io/ez_setup.py
sudo /usr/local/bin/python2.7 ez_setup.py
sudo /usr/local/bin/easy_install-2.7 pip

