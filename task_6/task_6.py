import numpy as np
from tools import trapeze_method, simpson_1_3_combined, simpson_3_8_combined, romberg_loop, romberg_recursion, open_quad

def function_1(x):
    return x * np.e ** (2 * x)

def function(x):
    return np.e ** (-(x ** 2))

limits = [0, 2]

dividing_steps = 20

analytic_result = 0.88208139

limits_1 = [0, 4]

analytic_result_1 = 5216.9265

if __name__ == '__main__':
    print(f'The analytic result is {analytic_result}\n')
    print(f'Using Trapeze method with {dividing_steps} steps, the result is {round(trapeze_method(function, limits, dividing_steps), 6)}\n')
    print(f'The difference between the analytic result and Trapeze method is {round(abs(analytic_result - trapeze_method(function, limits, dividing_steps)), 6)}\n')
    print(f'Using Simpson 1/3 combined method with {dividing_steps} steps, the result is {round(simpson_1_3_combined(function, limits, dividing_steps), 6)}\n')
    print(f'The difference between the analytic result and Simpson 1/3 method is {round(abs(analytic_result - simpson_1_3_combined(function, limits, dividing_steps)), 6)}\n')
    result, iteration = simpson_3_8_combined(function_1, limits_1, analytic_result_1)
    print(f'Using Simpson 1/3 combined method, the result is {round(result, 6)} after {iteration} iterations\n')
    print(f'My computer could not finish the calculations with the desired error (10 ** -5) so the best I managed to do is (10 ** -3). \nI wrote two Romberg functions to test it:\n')
    result_1, iteration_1 = romberg_loop(function_1, limits_1, analytic_result_1)
    print(f'The result using Romberg structured with a loop is {round(result_1, 6)} with {iteration_1} loops\n')
    result_2, iteration_2 = romberg_recursion(function_1, limits_1, analytic_result_1)
    print(f'The result using recursion Romberg is {round(result_2, 6)} with {iteration_2} iterations\n')
    result_3, iteration_3 = open_quad(function_1, limits_1, analytic_result_1)
    print(f'The result using open Quad method is {round(result_3, 6)} with {iteration_3} iterations\n')

