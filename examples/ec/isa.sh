#!/bin/sh
#删除已有文件

# for datablock in {10..10}
# do
#   for parityblock in {5..6}
#   do
#       rm isa_object_${datablock}_${parityblock}_1_avx2.txt
#   done
# done

# for datablock in {4..10}
# do
#   for parityblock in {2..4}
#   do
#     for((e=1;e<$parityblock;e++))
#     do 
#       rm isa_${datablock}_${parityblock}_${e}_1.txt
#     done
#     #e=${parityblock}
#   done
# done
# #开始保存运行结果
# for datablock in {4..10}
# do
#   for parityblock in {2..4}
#   do
#     for((e=1;e<$parityblock;e++))
#     do 
#       for i in {1..1}
#       do
#           ./ec_simple_example  -k $datablock -p $parityblock -e $e >> isa_${datablock}_${parityblock}_${e}_1.txt
#           sleep 10
#       done
#     done
#   done
# done

# #开始保存运行结果
# for datablock in {10..10}
# do
#   for parityblock in {5..6}
#   do
#     for i in {1..10}
#     do
#         ./ec_simple_example  -k $datablock -p $parityblock  >> isa_object_${datablock}_${parityblock}_1_avx2.txt
#         sleep 1
#     done
#   done
# done

#!/bin/sh
#删除已有文件


# #开始保存运行结果
# datablock=4
# parityblock=2 
# for i in {1..10}
# do
#     ./ec_simple_example  -k $datablock -p $parityblock  >> isa_object_${datablock}_${parityblock}_1_avx2.txt
#     sleep 1
# done
# datablock=6
# parityblock=3 
# for i in {1..10}
# do
#     ./ec_simple_example  -k $datablock -p $parityblock  >> isa_object_${datablock}_${parityblock}_1_avx2.txt
#     sleep 1
# done

# datablock=10
# parityblock=4 
# for i in {1..10}
# do
#     ./ec_simple_example  -k $datablock -p $parityblock  >> isa_object_${datablock}_${parityblock}_1_avx2.txt
#     sleep 1
# done


for datablock in {4..10}
do
  for parityblock in {2..4}
  do
    e=$parityblock
    rm ./data/isa_decode_${datablock}_${parityblock}_${e}.txt
  done
done
#开始保存运行结果
for datablock in {4..10}
do
echo $datablock
  for parityblock in {2..4}
  do
    for i in {1..1}
    do
      e=$parityblock
      ./ec_simple_example  -k $datablock -p $parityblock -e $e >> ./data/isa_decode_${datablock}_${parityblock}_${e}.txt
      sleep 10
    done
  done
done


