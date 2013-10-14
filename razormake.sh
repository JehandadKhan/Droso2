cd ./Release
module unload intel
module load intel/13.1.0
source /share/apps/intel/composer_xe_2013.2.146/composer_xe_2013.2.146/bin/compilervars.sh intel64
make -k all
cd ..
