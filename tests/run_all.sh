Test #1: 1D Dam=break with a triangular barrier
echo ====================================================================
echo
echo RUNNING TEST1: 1D DAMBREAK
echo
echo ====================================================================
echo
echo
cd test1_db1d
rm -r data results vis
rm xsurf
make
python setrun.py
export INDIR=data
./xsurf
python ver_1d.py
cd ..

echo ====================================================================
echo
echo RUNNING TEST2: 2D PARTIAL DAMBREAK
echo
echo ====================================================================
echo
echo
cd test2_db2d
rm -r data results vis
rm xsurf
make
python setrun.py
export INDIR=data
./xsurf
python ver.py
cd ..


echo ====================================================================
echo
echo RUNNING TEST3: Solitary wave over horizontal channel
echo
echo ====================================================================
echo
echo
cd test3_solitarywave
rm -r data results vis
rm xsurf
rm *.o
make
python setrun.py
export INDIR=data
./xsurf
python ver_1d.py
cd ..

echo ====================================================================
echo
echo RUNNING TEST4: Solitary wave over a sloped channel
echo
echo ====================================================================
echo
echo
cd test4_solitarysloped
rm -r data results vis
rm xsurf
rm *.o
make
python setrun.py
export INDIR=data
./xsurf
python ver_1d.py
cd ..

echo ====================================================================
echo
echo RUNNING TEST 5: 
echo Tsunami runup onto a complex three-dimensional beach Monai Valley
echo
echo ====================================================================
echo
echo
cd test5_labokushiri
rm -r data results vis
rm xsurf
rm *.o
make
python setrun.py
export INDIR=data
./xsurf
python ver.py
python gauges.py
cd ..