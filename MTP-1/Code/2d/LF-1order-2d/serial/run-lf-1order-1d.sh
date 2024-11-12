rm *.txt
rm *.x
rm *.exe
rm *.out

g++ lf-1order-2d.cpp -o finite-difference-grid-lf.x
./finite-difference-grid-lf.x

python plot-lf-1order-2d.py