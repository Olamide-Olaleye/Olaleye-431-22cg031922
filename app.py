from flask import Flask, render_template, request
import math

app = Flask(__name__)

# --- Helper Functions ---

def eval_poly(coeffs, x):
    """Evaluates polynomial at x using Horner's Method."""
    result = 0.0
    for c in coeffs:
        result = result * x + c
    return result

def get_poly_derivative(coeffs):
    """Computes derivative coefficients automatically."""
    n = len(coeffs) - 1
    if n < 1: return [0.0]
    deriv = []
    for i in range(n):
        power = n - i
        deriv.append(coeffs[i] * power)
    return deriv

def parse_coeffs(coeffs_str):
    try:
        return [float(x) for x in coeffs_str.strip().split()]
    except:
        return None

# --- Method Logic (Returns list of dicts for the table) ---

def solve_method(method, coeffs, x0, x1, delta, tol, max_iter):
    iterations = []
    root = None
    converged = False
    error_msg = None

    # Setup based on method
    curr_x = x0
    prev_x = x0 # For methods requiring history
    curr_b = x1 # For bracketing (upper bound)

    # Pre-calculate derivative for Newton
    deriv_coeffs = get_poly_derivative(coeffs) if method == 'newton' else []

    for i in range(1, max_iter + 1):
        row = {'iter': i}
        
        try:
            if method == 'bisection':
                # x0 is a, x1 is b
                fa = eval_poly(coeffs, x0)
                fb = eval_poly(coeffs, x1)
                
                if fa * fb >= 0:
                    error_msg = "Root not bracketed. f(a) and f(b) must have opposite signs."
                    break
                
                x_mid = (x0 + x1) / 2
                fx_mid = eval_poly(coeffs, x_mid)
                error = abs(x1 - x0)
                
                row.update({'a': x0, 'b': x1, 'x_val': x_mid, 'f_x': fx_mid, 'error': error})
                
                if abs(fx_mid) < tol or error < tol:
                    root = x_mid
                    converged = True
                    iterations.append(row)
                    break
                
                if fa * fx_mid < 0:
                    x1 = x_mid
                else:
                    x0 = x_mid
                    
            elif method == 'regula_falsi':
                fa = eval_poly(coeffs, x0)
                fb = eval_poly(coeffs, x1)
                
                if fa * fb >= 0:
                    error_msg = "Root not bracketed. f(a) and f(b) must have opposite signs."
                    break
                
                x_new = (x0 * fb - x1 * fa) / (fb - fa)
                fx_new = eval_poly(coeffs, x_new)
                error = abs(x_new - prev_x) if i > 1 else abs(x1 - x0) # First iter error approx
                
                row.update({'a': x0, 'b': x1, 'x_val': x_new, 'f_x': fx_new, 'error': error})
                
                if abs(fx_new) < tol:
                    root = x_new
                    converged = True
                    iterations.append(row)
                    break
                
                if fa * fx_new < 0:
                    x1 = x_new
                else:
                    x0 = x_new
                prev_x = x_new

            elif method == 'secant':
                # x0 is prev, x1 is curr
                f0 = eval_poly(coeffs, x0)
                f1 = eval_poly(coeffs, x1)
                
                if (f1 - f0) == 0:
                    error_msg = "Division by zero (f(x1) - f(x0) = 0)."
                    break
                
                x2 = x1 - (f1 * (x1 - x0)) / (f1 - f0)
                error = abs(x2 - x1)
                f2 = eval_poly(coeffs, x2)
                
                row.update({'x_prev': x0, 'x_curr': x1, 'x_val': x2, 'f_x': f2, 'error': error})
                
                if error < tol:
                    root = x2
                    converged = True
                    iterations.append(row)
                    break
                
                x0 = x1
                x1 = x2

            elif method == 'newton':
                fx = eval_poly(coeffs, curr_x)
                dfx = eval_poly(deriv_coeffs, curr_x)
                
                if dfx == 0:
                    error_msg = "Derivative is zero."
                    break
                
                x_next = curr_x - (fx / dfx)
                error = abs(x_next - curr_x)
                
                row.update({'x_curr': curr_x, 'f_x': fx, 'df_x': dfx, 'x_val': x_next, 'error': error})
                
                if error < tol:
                    root = x_next
                    converged = True
                    iterations.append(row)
                    break
                
                curr_x = x_next

            elif method == 'fixed_point':
                # Here coeffs represent g(x)
                x_next = eval_poly(coeffs, curr_x)
                error = abs(x_next - curr_x)
                
                row.update({'x_curr': curr_x, 'x_val': x_next, 'error': error})
                
                if error < tol:
                    root = x_next
                    converged = True
                    iterations.append(row)
                    break
                
                curr_x = x_next

            elif method == 'modified_secant':
                fx = eval_poly(coeffs, curr_x)
                fx_delta = eval_poly(coeffs, curr_x + delta * curr_x)
                
                denom = fx_delta - fx
                if denom == 0:
                    error_msg = "Division by zero in Mod Secant."
                    break
                
                x_next = curr_x - ((delta * curr_x * fx) / denom)
                error = abs(x_next - curr_x)
                
                row.update({'x_curr': curr_x, 'f_x': fx, 'x_val': x_next, 'error': error})
                
                if error < tol:
                    root = x_next
                    converged = True
                    iterations.append(row)
                    break
                
                curr_x = x_next

            iterations.append(row)

        except Exception as e:
            error_msg = f"Calculation Error: {str(e)}"
            break

    if not converged and not error_msg:
        root = iterations[-1]['x_val']
        error_msg = "Max iterations reached."

    return iterations, root, converged, error_msg

# --- Routes ---

@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    inputs = {}
    
    if request.method == 'POST':
        try:
            # Get basic inputs
            inputs = {
                'method': request.form.get('method'),
                'coeffs_str': request.form.get('coeffs'),
                'x0': float(request.form.get('x0', 0)),
                'x1': float(request.form.get('x1', 0) or 0), # Handle empty string
                'delta': float(request.form.get('delta', 0.01) or 0.01),
                'tol': float(request.form.get('tol', 0.0001)),
                'max_iter': int(request.form.get('max_iter', 10))
            }
            
            coeffs = parse_coeffs(inputs['coeffs_str'])
            if not coeffs:
                return render_template('index.html', error="Invalid coefficients.", inputs=inputs)

            iterations, root, converged, msg = solve_method(
                inputs['method'], coeffs, inputs['x0'], inputs['x1'], 
                inputs['delta'], inputs['tol'], inputs['max_iter']
            )
            
            result = {
                'iterations': iterations,
                'root': root,
                'converged': converged,
                'msg': msg
            }

        except ValueError:
            return render_template('index.html', error="Invalid numerical input.", inputs=inputs)

    return render_template('index.html', result=result, inputs=inputs)

if __name__ == '__main__':
    app.run(debug=True)