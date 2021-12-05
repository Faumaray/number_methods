use methods::{simple_gaus, EPS};
use std::{io, ops::Add, vec};
fn main() {
    let simple_iter_answers = simple_iter();
    let newton_answers = newton();
    for i in 0..simple_iter_answers.len() {
        println!(
            "Norma={:.5}\t Iteration={}",
            &simple_iter_answers[i].2, &simple_iter_answers[i].1
        );
        for j in 0..simple_iter_answers[i].0.len() {
            print!("x{}={:.5}\t", j, &simple_iter_answers[i].0[j]);
        }
        println!("");
    }
    println!("--------------------------------------------------------------------------------");
    println!("Newton:");
    for i in 0..newton_answers.len() {
        println!("Iteration={}", &newton_answers[i].1);
        for j in 0..newton_answers[i].0.len() {
            print!("x{}={:.5}\t", j, &newton_answers[i].0[j]);
        }
        println!("");
    }
}
fn norma(matrix: &Vec<Vec<f64>>, n: usize) -> f64 {
    let mut res = 0.0;
    for i in 0..n {
        for j in 0..n {
            res += matrix[i][j].powi(2);
        }
    }
    res = res.sqrt();
    return res;
}
fn func_iter(x: &Vec<f64>, i: usize) -> f64 {
    let mut res = 0.0;
    match i {
        0 => {
            res = 1.0 + ((x[1] - 0.5).sin() / 2.0);
        }
        1 => {
            res = 1.2 - (x[0].cos());
        }
        _ => {}
    }
    res
}
fn jacobi_matrix_iter(x: &Vec<f64>, i: usize, j: usize) -> f64 {
    let mut res = 0.0;
    match i {
        0 => match j {
            0 => {
                res = x[0].sin();
            }
            1 => {
                res = 0.0;
            }
            _ => {}
        },
        1 => match j {
            0 => {
                res = 0.0;
            }
            1 => {
                res = ((x[1] - 0.5).cos()) / 2.0;
            }
            _ => {}
        },
        _ => {}
    }
    res
}
fn simple_iter() -> Vec<(Vec<f64>, usize, f64)> {
    let n: usize = 2;
    let mut iter: usize = 0;
    let mut x0 = vec![0.0; n];
    let mut x = vec![0.0; n];
    let mut answers: Vec<(Vec<f64>, usize, f64)> = Vec::new();
    let mut matrix = vec![vec![0.0; n]; n];
    let mut max = 0.0;
    loop {
        for i in 0..n {
            for j in 0..n {
                matrix[i][j] = jacobi_matrix_iter(&x0, i, j);
            }
        }
        for i in 0..n {
            x[i] = func_iter(&x, i);
        }
        max = (x[0] - x0[0]).abs();
        for i in 1..n {
            if (x[i] - x0[i]).abs() > max {
                max = (x[i] - x0[i]).abs();
            }
        }
        let norma = norma(&matrix, n);
        x0 = x.clone();
        iter += 1;
        answers.push((x0.clone(), iter, norma));
        if (max < EPS) || (iter > 20) {
            break;
        }
    }
    answers
}

fn func_newton(x: &Vec<f64>, i: usize) -> f64 {
    let mut res = 0.0;
    match i {
        0 => {
            res = x[0].cos() + x[1] - 1.2;
        }
        1 => {
            res = 2.0 * x[0] - (x[1] - 0.5).sin() - 2.0;
        }
        _ => {}
    }
    res
}
fn jacobi_matrix_newton(x: &Vec<f64>, i: usize, j: usize) -> f64 {
    let mut res = 0.0;
    match i {
        0 => match j {
            0 => {
                res = -(x[0].sin());
            }
            1 => {
                res = 1.0;
            }
            _ => {}
        },
        1 => match j {
            0 => {
                res = 2.0;
            }
            1 => {
                res = -((x[1] - 0.5).cos());
            }
            _ => {}
        },
        _ => {}
    }
    res
}
fn newton() -> Vec<(Vec<f64>, usize)> {
    let n: usize = 2;
    let mut iter: usize = 0;
    let mut dx = vec![0.0; n];
    let mut x = vec![0.0; n];
    let mut f = vec![0.0; n];
    let mut answers: Vec<(Vec<f64>, usize)> = Vec::new();
    let mut a = vec![vec![0.0; n]; n];
    let mut max: f64 = 0.0;
    loop {
        for i in 0..n {
            for j in 0..n {
                a[i][j] = jacobi_matrix_newton(&x, i, j);
            }
        }
        for i in 0..n {
            f[i] = -1.0 * func_newton(&x, i)
        }
        dx = simple_gaus(a.clone(), f.clone());
        max = dx[0].abs();
        for i in 1..n {
            if dx[i].abs() > max {
                max = dx[i].abs();
            }
        }
        for i in 0..n {
            x[i] += dx[i];
        }
        iter += 1;
        answers.push((x.clone(), iter));
        if max <= EPS || iter > 5 {
            break;
        }
    }
    answers
}
