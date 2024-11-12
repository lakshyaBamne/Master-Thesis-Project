rm *.txt
rm *.x
rm *.exe
rm *.out

g++ lw-2order-1d.cpp -o finite-difference-grid-lw.x
./finite-difference-grid-lw.x

python plot-lw-2order-1d.py