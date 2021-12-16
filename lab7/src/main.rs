use std::{io, ops::Add, vec};

fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut n = usize::MIN;
    println!("Ручной ввод?\n(1)Да;\n(2)Нет");
    let mut choice = String::new();
    io::stdin().read_line(&mut choice).expect("Failed to read");
    let mut inp = 0;
    match choice.trim().parse::<usize>() {
        Ok(v) => {
            inp = v;
        }
        Err(..) => eprintln!("Не подходящий формат: {}", choice.trim()),
    };
    if inp == 1 {
        println!("Введите размерность для квадратной матрицы");
        let mut input_text = String::new();
        io::stdin()
            .read_line(&mut input_text)
            .expect("Failed to read");
        let trimmed = input_text.trim();
        match trimmed.parse::<usize>() {
            Ok(v) => {
                n = v;
            }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
        };
        for _i in 0..n {
            let mut temp = Vec::<f64>::new();
            for _j in 0..n {
                temp.push(0.0);
            }
            matrix.push(temp);
        }
        println!("Матрицу: ");
        for i in 0..n {
            for j in 0..n {
                let mut input_text = String::new();
                io::stdin()
                    .read_line(&mut input_text)
                    .expect("Failed to read");
                let trimmed = input_text.trim();
                match trimmed.parse::<f64>() {
                    Ok(v) => {
                        matrix[i][j] = v;
                    }
                    Err(..) => println!("Не подходящий формат: {}", trimmed),
                };
            }
        }
    }
    let y = vec![
        0.17, 0.07, 0.17, 0.05, 0.12, 0.00, 0.01, -0.05, -0.21, -0.50, -0.50, -0.86, -1.24, -1.47,
        -1.79, -2.25, -2.55, -3.18, -3.60, -3.93,
    ];
    let x = vec![
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        17.0, 18.0, 19.0, 20.0,
    ];
    n = 20;
    println!("Kvadrat = {} for x={}", kvadr(n, &x, &y, 3.0), 3.0);
    println!("Newton = {} for x={}", newton(n, &x, &y, 20.0), 21.0);
    println!(
        "Lagr = {} for x={}",
        lagr(n, x.clone(), y.clone(), 20.0),
        21.0
    );
    println!("Spline = {} for x={}", spline(n, x, y, 20.0), 21.0);
}
fn lagr(n: usize, x: Vec<f64>, y: Vec<f64>, q: f64) -> f64 {
    let mut result = 0.0;
    let mut s = 0.0;
    for i in 0..n {
        s = 1.0;
        for j in 0..n {
            if j != i {
                s *= (q - x[j]) / (x[i] - x[j]);
            }
        }
        result += y[i] * s;
    }
    result
}
#[derive(Default, Clone, Copy)]
struct SplineTuple {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    x: f64,
}
fn spline(n: usize, x: Vec<f64>, y: Vec<f64>, q: f64) -> f64 {
    let mut s = SplineTuple::default();
    let splines = build_spline(n, &x, &y);
    if q <= splines[0].x {
        s = splines[0];
    } else if q >= splines[n - 1].x {
        s = splines[n - 1];
    } else {
        let mut i: usize = 0;
        let mut j: usize = n - 1;
        while i + 1 < j {
            let mut k = i + (j - i) / 2;
            if q <= splines[k].x {
                j = k;
            } else {
                i = k;
            }
        }
        s = splines[j];
    }
    let dx = q - s.x;
    (s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx)
}
fn build_spline(n: usize, x: &Vec<f64>, y: &Vec<f64>) -> Vec<SplineTuple> {
    let mut splines = vec![SplineTuple::default(); n];
    for i in 0..n {
        splines[i].x = x[i];
        splines[i].a = y[i];
    }
    splines[0].c = 0.0;
    splines[n - 1].c = 0.0;
    let mut alpha = vec![0.0; n - 1];
    let mut beta = vec![0.0; n - 1];
    alpha[0] = 0.0;
    beta[0] = 0.0;
    for i in 1..n - 1 {
        let hi = x[i] - x[i - 1];
        let hi_one = x[i + 1] - x[i];
        let a = hi;
        let c = 2.0 * (hi + hi_one);
        let b = hi_one;
        let f = 6.0 * ((y[i + 1] - y[i]) / hi_one - (y[i] - y[i - 1]) / hi);
        let z = (a * alpha[i - 1] + c);
        alpha[i] = -b / z;
        beta[i] = (f - a * beta[i - 1]) / z;
    }

    for i in (0..(n - 2)).rev() {
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
    }

    for i in (1..=(n - 1)).rev() {
        let hi = x[i] - x[i - 1];
        splines[i].d = (splines[i].c - splines[i - 1].c) / hi;
        splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi;
    }
    splines
}
fn newton(n: usize, x: &Vec<f64>, y: &Vec<f64>, q: f64) -> f64 {
    let mut result = 0.0;
    let mut m = vec![vec![0.0; n * 2]; n * 2];
    for i in 1..n {
        m[0][i - 1] = y[i] - y[i - 1];
    }
    for i in 1..n {
        let mut f = 1.0;
        for k in (1..=(i + 1)).rev() {
            f *= k as f64;
        }
        for j in (1..(n - 1)).rev() {
            m[i][j - 1] = (m[i - 1][j] - m[i - 1][j - 1]) / f;
        }
    }
    if x[n - 1] - q - q - x[1] < 0.0 {
        let mut s = y[n - 1];
        for i in 0..n - 1 {
            let mut p = 1.0;
            for j in 0..=i {
                p *= q - ((n - j) as f64);
            }
            s += p * m[i][n - i - 2];
        }
        result = s;
    } else {
        let mut s = y[0];
        for i in 0..n - 1 {
            let mut p = 1.0;
            for j in 0..=i {
                p *= q - x[j];
            }
            s += p * m[i][0];
        }
        result = s;
    }

    result
}
fn kvadr(n: usize, x: &Vec<f64>, y: &Vec<f64>, q: f64) -> f64 {
    let mut result = 0.0;
    let mut a = vec![vec![0.0; 3]; 3];
    let mut b = vec![0.0; 3];
    let mut c = [0.0; 3];
    a[0][0] = n as f64;
    for i in 0..n {
        for j in 0..3 {
            for k in 0..3 {
                let power = j + k;
                a[j][k] += x[i].powi(power as i32);
            }
            b[j] += y[i] * x[i].powi(j as i32);
        }
    }
    b = methods::simple_gaus(a, b);
    result = b[0] + (b[1] * q) + (b[2] * q * q);
    result
}
