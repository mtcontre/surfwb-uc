for TEST in test1_db1d test2_db2d test3_solitarywave test4_solitarysloped test5_labokushiri
do
  cd $TEST 
#   eog vis
  rm -r data results vis
  rm xsurf
  rm *.o
  cd ..
done


