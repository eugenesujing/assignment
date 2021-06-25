Two versions have been provided, **source_int.cpp** takes *int* as input while **source_double.cpp** takes *double* as input.

For **source_int.cpp**:
1. Compile it using the following command in your terminal:
`g++ -std=c++11 source_int.cpp -o int -Wall`
2. Run the programm by using `./int [G value] [B value] [R value]`. For example:
`./int 255 255 255`
The output would be like this :
```
It took 2.52644s for direct multiplication
It took 1.79667s for lift_based calculation
result of direct_conversion():
Cg = 0  Y = -3 Co = 0
result of lift_conversion():
Cg = 0  Y = -2 Co = 0
Different result for channel Y
```
There will be notification if the results from two methods are inconsistent. And bad input will terminate the programm immediately.


For **source_double.cpp**:
1. Compile it using the following command in your terminal:
`g++ -std=c++11 source_double.cpp -o double -Wall`
2. Run the programm by using `./double [G value] [B value] [R value]`. For example:
`./int 255.0 255.0 255.0`
The output would be like this:
```
It took 0.0341234s for direct multiplication
It took 0.0211219s for lift_based calculation
result of direct_conversion():
Cg = 57.5  Y = 197.5 Co = -115
result of lift_conversion():
Cg = 57.5  Y = 197.5 Co = -115

```
