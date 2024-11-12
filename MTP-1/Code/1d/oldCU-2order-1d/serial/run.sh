rm *.txt
rm *.x
rm *.exe
rm *.out

g++ main.cpp -o oldCU.x
./oldCU.x

python plot.py