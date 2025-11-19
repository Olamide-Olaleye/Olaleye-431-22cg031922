import sys

# --- Helpers ---
def get_input(prompt, type_func=float):
    while True:
        try: return type_func(input(prompt).strip())
        except ValueError: print("Invalid input.")

def get_coeffs(prompt):
    while True:
        try: return [float(x) for x in input(prompt).split()]
        except ValueError: print("Invalid coeffs.")

def eval_poly(c, x):
    r = 0.0
    for a in c: r = r * x + a
    return r

def poly_str(c):
    n = len(c) - 1
    return " + ".join(f"{v:.2f}x^{n-i}" if n-i > 0 else f"{v:.2f}" for i, v in enumerate(c) if v != 0) or "0"

def solve_loop(headers, state, step_func, tol, max_iter, coeffs=None):
    print(f"\n{' | '.join(f'{h:^12}' for h in headers)}")
    print("-" * (15 * len(headers)))
    
    for i in range(1, max_iter + 1):
        try:
            new_state, row, err, root = step_func(state)
        except ZeroDivisionError:
            print("Error: Division by zero.")
            return
        
        print(f"{i:^4} | " + " | ".join(f"{v:^12.6f}" for v in row))
        
        # Check convergence on error OR function value (if coeffs provided)
        f_val = abs(eval_poly(coeffs, root)) if coeffs else float('inf')
        if err < tol or f_val < tol:
            print(f"\nConverged: {root:.6f} in {i} iters.")
            return
        state = new_state
    print(f"\nMax iterations reached. Approx root: {root:.6f}")

# --- Method Generators ---
def bracketing_method(name, calc_x):
    print(f"\n--- {name} ---")
    c = get_coeffs("Coeffs (an ... a0): ")
    print(f"P(x) = {poly_str(c)}")
    a, b = get_input("Lower (a): "), get_input("Upper (b): ")
    
    if eval_poly(c, a) * eval_poly(c, b) >= 0: return print("Error: f(a)*f(b) must be < 0")
    
    def step(s):
        a, b = s['a'], s['b']
        fa, fb = eval_poly(c, a), eval_poly(c, b)
        xm = calc_x(a, b, fa, fb)
        fxm = eval_poly(c, xm)
        # If signs change between a and xm, root is there; else it's between xm and b
        new_s = {'a': a, 'b': xm} if fa * fxm < 0 else {'a': xm, 'b': b}
        return (new_s, [a, b, xm, fxm, abs(b-a)], abs(b-a), xm)

    solve_loop(["a", "b", "x_new", "f(x)", "Error"], {'a': a, 'b': b}, step, get_input("Tol: "), int(get_input("Iter: ")), c)

def open_method(name, step_logic, setup_fn=None):
    print(f"\n--- {name} ---")
    c = get_coeffs("Coeffs: ")
    print(f"Function: {poly_str(c)}")
    extra = setup_fn(c) if setup_fn else {} # Store derivative coeffs or delta
    x0 = get_input("Initial (x0): ")
    
    def step(s):
        x = s['x']
        xn = step_logic(c, x, extra)
        err = abs(xn - x)
        return ({'x': xn}, [x, xn, err], err, xn)

    solve_loop(["x_curr", "x_next", "Error"], {'x': x0}, step, get_input("Tol: "), int(get_input("Iter: ")), c)

def secant_wrapper():
    print("\n--- Secant Method ---")
    c = get_coeffs("Coeffs: ")
    x0, x1 = get_input("Guess 1 (x0): "), get_input("Guess 2 (x1): ")
    
    def step(s):
        x0, x1 = s['x0'], s['x1']
        f0, f1 = eval_poly(c, x0), eval_poly(c, x1)
        x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0)
        return ({'x0': x1, 'x1': x2}, [x0, x1, x2, abs(x2-x1)], abs(x2-x1), x2)

    solve_loop(["x_prev", "x_curr", "x_next", "Error"], {'x0': x0, 'x1': x1}, step, get_input("Tol: "), int(get_input("Iter: ")), c)

# --- Main ---
def main():
    # Map logic to readable keys
    methods = {
        '1': lambda: bracketing_method("Bisection", lambda a,b,fa,fb: (a+b)/2),
        '2': lambda: bracketing_method("Regula Falsi", lambda a,b,fa,fb: (a*fb - b*fa)/(fb-fa)),
        '3': secant_wrapper,
        '4': lambda: open_method("Newton-Raphson", 
                                 lambda c,x,d: x - eval_poly(c,x)/eval_poly(d,x), 
                                 # Lambda to compute derivative coefficients
                                 lambda c: [c[i]*(len(c)-1-i) for i in range(len(c)-1)]),
        '5': lambda: open_method("Fixed Point", lambda c,x,_: eval_poly(c,x)),
        '6': lambda: open_method("Modified Secant", 
                                 lambda c,x,d: x - (d*x*eval_poly(c,x))/(eval_poly(c,x+d*x)-eval_poly(c,x)), 
                                 lambda _: get_input("Delta: "))
    }
    
    while True:
        print("\n=== ZOF Poly Solver ===")
        print("1.Bisection  2.Regula Falsi  3.Secant\n4.Newton     5.Fixed Point   6.Mod Secant  Q.Quit")
        ch = input("Select: ").strip().upper()
        if ch == 'Q': sys.exit()
        if ch in methods: methods[ch]()
        else: print("Invalid selection.")

if __name__ == "__main__": main()