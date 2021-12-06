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
    } else {
        matrix = vec![
            vec![0.58993, 0.036004, -0.32946, 0.27315],
            vec![-4.3515, -0.24076, 2.2030, -2.0148],
            vec![6.7688, 0.39005, -3.5692, 3.1341],
            vec![-31.085, -1.7198, 15.737, -14.393],
        ];
        n = 4;
    }
    let power_ans = power_iteration(&matrix, 100);
    for value in power_ans {
        print!("|{:.06}|\t", value);
    }
    println!("\nQR");
    let ans = qr(matrix);
    for row in ans {
        for value in row {
            print!("|{:.06}|\t", value);
        }
        println!("");
    }
}

fn qr(matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let mut out = vec![vec![0.0; n]; n];
    let mut q = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                q[i][j] = 1.0;
            } else {
                q[i][j] = 0.0;
            }
        }
    }
    let mut qh = vec![vec![0.0; n]; n];
    let mut r = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            r[i][j] = matrix[i][j];
        }
    }
    let mut qr_a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            qr_a[i][j] = matrix[i][j];
        }
    }
    let mut qrk: usize = 0;
    let mut k: usize = 0;
    while qrk < 10 {
        for i in 0..n {
            for j in 0..n {
                r[i][j] = qr_a[i][j];
            }
        }
        k = 0;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    q[i][j] = 1.0;
                } else {
                    q[i][j] = 0.0;
                }
            }
        }

        while k < (n - 1) {
            let mut a = vec![0.0; n - k];
            for i in k..n {
                a[i - k] = r[i][k];
            }
            let u_vec = reflection_for_vector(&a, n - k);
            let u_matrix = reflection_for_matrix(&u_vec, n - k);

            qh = add_matrix(&u_matrix, k, n);

            q = matrix_multi(&q, &qh);

            r = matrix_multi(&qh, &r);

            k += 1;
        }

        qr_a = matrix_multi(&r, &q);

        for i in 0..n {
            for j in 0..n {
                out[i][j] = qr_a[i][j];
            }
        }
        qrk += 1;
    }
    out
}
fn vector_normir(a: &Vec<f64>, n: usize) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        sum += a[i].powi(2);
    }
    sum.sqrt()
}
fn reflection_for_vector(a: &Vec<f64>, n: usize) -> Vec<f64> {
    let mut u = vec![0.0; n];
    u[0] = a[0] - vector_normir(&a, n);
    for i in 1..n {
        u[i] = a[i];
    }
    u
}
fn reflection_for_matrix(u: &Vec<f64>, n: usize) -> Vec<Vec<f64>> {
    let mut refl_matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                refl_matrix[i][j] = 1.0 - 2.0 * u[i] * u[j] / vector_normir(&u, n).powi(2);
            } else {
                refl_matrix[i][j] = -2.0 * u[i] * u[j] / vector_normir(&u, n).powi(2);
            }
        }
    }
    refl_matrix
}
fn add_matrix(u: &Vec<Vec<f64>>, k: usize, n: usize) -> Vec<Vec<f64>> {
    let mut full_matrix = vec![vec![0.0; n]; n];
    for i in 0..k {
        for j in 0..k {
            if i == j {
                full_matrix[i][j] = 1.0;
            } else {
                full_matrix[i][j] = 0.0;
            }
        }
    }
    for i in k..n {
        for j in k..n {
            full_matrix[i][j] = u[i - k][j - k];
        }
    }
    full_matrix
}
fn matrix_multi(first: &Vec<Vec<f64>>, second: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = first.len();
    let mut result = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                result[i][j] += first[i][k] * second[k][j];
            }
        }
    }
    result
}
fn matr_umn_mas(matrix: &Vec<Vec<f64>>, b: &Vec<f64>, n: usize, m: usize, k: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    for i in 0..n {
        for j in 0..k {
            let mut s = 0.0;
            for l in 0..m {
                s += matrix[i][l] * b[l];
            }
            out[i] = s;
        }
    }
    out
}
fn norm(vec: &Vec<f64>) -> f64 {
    let mut s = 0.0;
    for i in 0..vec.len() {
        s += vec[i].powi(2);
    }
    s.sqrt()
}
fn nib(vc: &Vec<f64>, s: f64) -> Vec<f64> {
    let mut out = vec![0.0; vc.len()];
    for i in 0..vc.len() {
        out[i] = vc[i] / s;
    }
    out
}
use rand::prelude::*;
fn power_iteration(matrix: &Vec<Vec<f64>>, iter: usize) -> Vec<f64> {
    let mut rng = thread_rng();
    let mut b_k: Vec<f64> = vec![rng.gen(); matrix.len() * matrix.len()];
    for _i in 0..iter {
        let mut b_k1 = matr_umn_mas(matrix, &b_k, matrix.len(), matrix.len(), matrix.len());
        let mut b_k1_norm = norm(&b_k1);
        b_k = nib(&b_k1, b_k1_norm);
    }

    b_k
}
