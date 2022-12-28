#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <stdlib.h>

#define ESC 27

// create the structure that fully describe complexequation
typedef struct complex_equations
{
    double real_in_numerator;
    double real_in_denominator;
    double imeginary_in_numerator;
    double imeginary_in_denominator;
    double real_in_answer;
    double imeginary_in_answer;
} complex_equations;

int verify_int(int *);
int verify_choise_case(int);
int verify_double(double *);
int verify_bounds(double, double, double *, double *);
int verify_double_greater_then_zero(double);
int verify_step(double, double);
void clear();
void assign_variables_for_the_first_complex_equation(double, double, double, double, double);
void assign_variables_for_the_second_complex_equation(double, double, double, double, double);
void assign_variables_for_the_third_complex_equation(double, double, double, double, double);
void assign_variables_for_the_fourth_complex_equation(double, double, double, double, double);
void calculate_complex_resistance(complex_equations);
void compute_the_table(double, double, double, double, double, double, int, double, void (*)(double, double, double, double, double), complex_equations);
void display_anwer_for_complex_equation(complex_equations);

int main()
{
    char key;
    do
    {
        double lower_bound, upper_bound, first_bound, second_bound, result, resistance_1, resistance_2, inductance, electric_capacity, step, omega;
        int choise_scheme, e;

        system("cls");

        // Start entering data! >>>>
        do
        {
            printf("Enter the scheme number: ");
        } while ((verify_int(&choise_scheme) == 1) || (verify_choise_case(choise_scheme) == 1));

        // we have at least one resistor
        do
        {
            printf("Enter thr resistance of the first resistor, R1: ");
        } while ((verify_double(&resistance_1) == 1) || (verify_double_greater_then_zero(resistance_1) == 1));

        // in case 3th and 4th scheme we have second resistor /start_1
        if (choise_scheme > 2)
        {
            do
            {
                printf("Enter thr resistance of the second resistor, R2: ");
            } while ((verify_double(&resistance_2) == 1) || (verify_double_greater_then_zero(resistance_1) == 1));
        }
        // finish_1

        // keep entering data and verify that they greater then zero
        do
        {
            printf("Enter inductance, L: ");
        } while ((verify_double(&inductance) == 1) || (verify_double_greater_then_zero(inductance) == 1));
        do
        {
            printf("Enter electric capacity, C: ");
        } while ((verify_double(&electric_capacity) == 1) || (verify_double_greater_then_zero(electric_capacity) == 1));

        // we enter the border in arbitrarily way, function verify_bounds decide whic one is greater /start_2
        do
        {
            printf("Enter the first bound: ");
        } while ((verify_double(&first_bound) == 1));
        do
        {
            printf("Enter the second bound: ");
        } while ((verify_double(&second_bound) == 1) || (verify_bounds(first_bound, second_bound, &lower_bound, &upper_bound)));
        // finish_2

        do
        {
            printf("Enter step: ");
        } while ((verify_double(&step) == 1) || (verify_step(step, upper_bound - lower_bound) == 1));
        // >>>>Stop entering data!

        switch (choise_scheme)
        {
        case 1:
            complex_equations first_equation;
            compute_the_table(electric_capacity, inductance, lower_bound, upper_bound, resistance_1, resistance_2, choise_scheme, step, assign_variables_for_the_first_complex_equation, first_equation);
            break;
        case 2:
            complex_equations second_equation;
            compute_the_table(electric_capacity, inductance, lower_bound, upper_bound, resistance_1, resistance_2, choise_scheme, step, assign_variables_for_the_second_complex_equation, second_equation);
            break;
        case 3:
            complex_equations third_equation;
            compute_the_table(electric_capacity, inductance, lower_bound, upper_bound, resistance_1, resistance_2, choise_scheme, step, assign_variables_for_the_third_complex_equation, third_equation);
            break;
        case 4:
            complex_equations fourth_equation;
            compute_the_table(electric_capacity, inductance, lower_bound, upper_bound, resistance_1, resistance_2, choise_scheme, step, assign_variables_for_the_fourth_complex_equation, fourth_equation);
            break;
        }

        // Ask to keep app work
        printf("Press any button to keep going or escape to quit.");
        key = getch();

    } while (key != ESC);
}

// idea is to have function that can assign variables seperatly for each complex equation
// date stored as a global varieble in stucture
// start of the block of assignment-------------->

// add resistance_2 to unify the look of functions that assign variables
void assign_variables_for_the_first_complex_equation(double inductance, double electric_capacity, double resistance_1, double resistance_2, double omega)
{
    complex_equations first_equation;
    // assign variables using given formulas
    first_equation.real_in_numerator = inductance / electric_capacity;
    first_equation.real_in_denominator = resistance_1;
    first_equation.imeginary_in_numerator = ((0 - resistance_1) / (omega * electric_capacity));
    first_equation.imeginary_in_denominator = (omega * inductance - 1 / (omega * electric_capacity));

    // pass the number/ name of structure of equation
    calculate_complex_resistance(first_equation);
}

void assign_variables_for_the_second_complex_equation(double inductance, double electric_capacity, double resistance_1, double resistance_2, double omega)
{
    complex_equations second_equation;
    // assign variables using given formulas
    second_equation.real_in_numerator = inductance / electric_capacity;
    second_equation.real_in_denominator = resistance_1;
    second_equation.imeginary_in_numerator = resistance_1 / (omega * electric_capacity);
    second_equation.imeginary_in_denominator = (omega * inductance - 1 / (omega * electric_capacity));

    // pass the number/ name of structure of equation
    calculate_complex_resistance(second_equation);
}

void assign_variables_for_the_third_complex_equation(double inductance, double electric_capacity, double resistance_1, double resistance_2, double omega)
{
    complex_equations third_equation;
    // assign variables using given formulas
    third_equation.real_in_numerator = resistance_1 * resistance_2;
    third_equation.real_in_denominator = resistance_1 * (omega * inductance - 1 / (omega * electric_capacity));
    third_equation.imeginary_in_numerator = resistance_1 + resistance_2;
    third_equation.imeginary_in_denominator = (omega * inductance - 1 / (omega * electric_capacity));

    // pass the number/ name of structure of equation
    calculate_complex_resistance(third_equation);
}

void assign_variables_for_the_fourth_complex_equation(double inductance, double electric_capacity, double resistance_1, double resistance_2, double omega)
{
    complex_equations fourth_equation;
    // assign variables using given formulas
    fourth_equation.real_in_numerator = resistance_1 * resistance_2 + inductance / electric_capacity;
    fourth_equation.real_in_denominator = resistance_1 * (omega * inductance * resistance_1 - resistance_2 / (omega * electric_capacity));
    fourth_equation.imeginary_in_numerator = resistance_1 + resistance_2;
    fourth_equation.imeginary_in_denominator = (omega * inductance - 1 / (omega * electric_capacity));

    // pass the number/ name of structure of equation
    calculate_complex_resistance(fourth_equation);
}
// finish of the block of assignment  ----------------<

//  function that can calculate the answer for specific name of structure
void calculate_complex_resistance(complex_equations number_equation)
{
    // using general formula calculate separately real number and imaginary one in answer
    number_equation.real_in_answer = (number_equation.real_in_numerator * number_equation.real_in_denominator + number_equation.imeginary_in_numerator * number_equation.imeginary_in_denominator) / (number_equation.real_in_denominator * number_equation.real_in_denominator + number_equation.imeginary_in_denominator * number_equation.imeginary_in_denominator);

    number_equation.real_in_answer = (number_equation.imeginary_in_numerator * number_equation.real_in_denominator - number_equation.real_in_numerator * number_equation.imeginary_in_denominator) / (number_equation.real_in_denominator * number_equation.real_in_denominator + number_equation.imeginary_in_denominator * number_equation.imeginary_in_denominator);
}

// idea is to calculate the table using loop that call function that calculate the specific name of structure
void compute_the_table(double electric_capacity, double inductance, double lower_bound, double upper_bound, double resistance_1, double resistance_2, int choise_scheme, double step, void (*scheme)(double, double, double, double, double), complex_equations number_equation)
{
    double resonance_frequency, omega;
    resonance_frequency = 1 / (2 * M_PI * sqrt(inductance * electric_capacity));
    printf("Resonance frequency: %e", resonance_frequency);
    printf("\nFrequency\tResistance\n");

    while (lower_bound < upper_bound)
    {
        omega = 2 * M_PI * lower_bound;
        (*scheme)(inductance, electric_capacity, resistance_1, resistance_2, omega);
        printf("%7g\t\t", lower_bound);
        display_anwer_for_complex_equation(number_equation);
        printf("\n");
        lower_bound += step;
    }
}

void display_anwer_for_complex_equation(complex_equations number_equation)
{
    if (number_equation.imeginary_in_answer >= 0)
    {
        printf("%g + i*%g", number_equation.real_in_answer, number_equation.imeginary_in_answer);
    }
    else
    {
        printf("%g - i*%g", number_equation.real_in_answer, fabs(number_equation.imeginary_in_answer));
    }
}

int verify_int(int *variable_int)
{
    char ch;
    scanf("%d%c", variable_int, &ch);
    if (ch == '\n')
    {
        return 0;
    }
    printf("Invalid number. Try again!\n");
    clear();
    return 1;
}

int verify_choise_case(int var_choise_case)
{
    if ((var_choise_case < 5) && (var_choise_case > 0))
    {
        return 0;
    }
    printf("Invalid choise of equation or method. Try again!\n");
    return 1;
}

int verify_double(double *variable_double)
{
    char ch;
    scanf("%lf%c", variable_double, &ch);
    if (ch == '\n')
    {
        return 0;
    }
    printf("Invalid number. Try again!\n");
    clear();
    return 1;
}

int verify_double_greater_then_zero(double num)
{
    if (num > 0)
    {
        return 0;
    }
    return 1;
}

int verify_bounds(double var_first, double var_second, double *var_low, double *var_top)
{
    if (var_first < var_second)
    {
        *var_low = var_first;
        *var_top = var_second;
        return 0;
    }
    else if (var_first > var_second)
    {
        *var_top = var_first;
        *var_low = var_second;
        return 0;
    }
    printf("Invalid borders. Try again!\n");
    return 1;
}

int verify_step(double step, double range)
{
    if (step > range)
    {
        return 1;
    }
    return 0;
}

void clear()
{
    int character;
    while ((character = getchar()) != '\n' && character != EOF) // end of file
        ;
}
