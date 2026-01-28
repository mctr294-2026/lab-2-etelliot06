#include <functional>

/* Tries to find a zero crossing in f() in the interval [a,b] with the bisection method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool bisection(std::function<double(double)> f, double a, double b, double *root) 
{
    // If there is a root, the functions will have opposite signs.
    if (f(a) * f(b) > 0)
        return false;

    for(int i=0; i < 100; ++i) // loop 100 times to make smaller intervals containing a root
    {
        double c = (a+b)/2; // find midpoint of the interval [a,b]
        if (std::abs(f(c)) < 1) // check if c is a root
        {
            *root = c;
            return true;
        }
        if (f(a) * f(c) < 0) // check if there is a root between [a,c]
        {
            b = c; // if true make the new interval [a,c]
        }
        else
        {
            a = c; // if false make the new interval [b,c]
        }
    }
    *root = (a+b)/2; // root is the located in the middle of the last interval
    return true;
}


/* Tries to find a zero crossing in f() in the interval [a,b] with the
 * false positive / regula falsi method
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found. Finding a root
 * is only guarenteed if f is continous within the interval and
 * a & b have opposite signs
 */
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root)
{
    // make a guess where a root can be using a straigh line from point A to B
    // same as bisection method we need the functions to have opposite signs
    if (f(a) * f(b) > 0)
        return false;
    
    double fa = f(a);
    double fb = f(b);

    for(int i=0; i < 100; ++i) // loop 100 times to make smaller intervals containing a root
    {
        double c = b - fb * (b-a) / (fb - fa); // the x value for the line
        double fc = f(c);

        if (std::abs(fc) < 1) // check if c is a root
        {
            *root = c;
            return true;
        }
        if (fa * fc < 0) // check if there is a root between [a,c]
        {
            b = c; // if true make the new interval [a,c]
        }
        else
        {
            a = c; // if false make the new interval [b,c]
        }
    }
    *root = b - fb * (b-a) / (fb - fa); // root is located at the line in the last interval
    return true;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the netwon-raphson method, given a function that computes the
 * derivative g() and a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or the derivative is zero.
 */
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root)
{
    for(int i=0; i < 100; ++i) // loop 100 times to make smaller intervals containing a root
    {
        double cnplusone = c - (f(c)/g(c)); // a tangent line is created at our guess, this formula shows where that tangent line is zero
        c = cnplusone;
        if (cnplusone < a || cnplusone > b) // if new value is outside the interval --> no root
            return false;
        if (abs(g(cnplusone)) < 0.001) // if the derivative approaches zero --> no root
            return false;
    }
    *root = c - (f(c)/g(c));
    return true;
}

/* Tries to find a zero crossing in f() in the interval [a,b] with
 * the secant method, given a starting guess c.
 * Returns true if a root is found. The crossing is stored in root.
 * Returns false if a crossing could not be found, which can happen
 * if iteration leaves the interval, or derivative is zero.
 */
bool secant(std::function<double(double)> f, double a, double b, double c, double *root)
{
    double cnminusone = 0.0; // the previous number has not been dertermined yet
    for(int i=0; i < 100; ++i) // loop 100 times to make smaller intervals containing a root
    {
        double cnplusone = c - (f(c)/((f(c)-f(cnminusone))/(c-cnminusone))); // a secant line is created at our guess, this formula shows where that secant line is zero
        cnminusone = c; // now that we ran it once, the new cnminusone becomes c
        c = cnplusone;
        if (cnplusone < a || cnplusone > b) // like newtons method, we look if its still in the interval
            return false;
        //if (abs(g(cnplusone)) < 0.001) // and if the derivative is still greater than zero BUT we don't have the derivative in this function
            //return false;
    }
    *root = c - (f(c)/((f(c)-f(cnminusone))/(c-cnminusone)));
    return true;
}