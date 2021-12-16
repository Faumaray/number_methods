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
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in &matrix {
        println!("{:?}", i);
    }
    println!("--------------------------------------------------------------------------------");
    println!("Laverre: ");
    println!("");
    let answer = laverre(matrix.clone());
    println!("Laverre eigenvals: ");
    for value in answer.0 {
        println!("{:?}", value);
    }
    println!("Laverre eigenvecs: ");
    for value in answer.1 {
        println!("{:?}", value);
    }
    println!("--------------------------------------------------------------------------------");
    println!("Krilov: ");
    println!("");
    let answer_krilov = krilov(matrix.clone());
    println!("Krilov eigenvals: ");
    for value in answer_krilov.0 {
        println!("{:?}", value);
    }
    println!("Krilov eigenvecs: ");
    for value in answer_krilov.1 {
        println!("{:?}", value);
    }
    println!("--------------------------------------------------------------------------------");
    println!("Fadeev: ");
    println!("");
    let answer_fadeev = fadeev(matrix);
    println!("Fadeev eigenvals: ");
    for value in answer_fadeev.0 {
        println!("{:?}", value);
    }
    println!("Fadeev eigenvecs: ");
    for value in answer_fadeev.1 {
        println!("{:?}", value);
    }
}

//проверяет входит ли в массив зна-чение Znach,
//используется при вычислении определителя
fn vkl(per: &Vec<usize>, znach: usize) -> bool {
    for value in per {
        if *value == znach {
            return true;
        }
    }
    false
}
//для определителя указывает знак с каким входит
//в сумму очередное слагаемое
fn perestanovka(per: &Vec<usize>) -> bool {
    let mut count = 0;
    let mut result = true;
    for i in 0..=(per.len() - 2) {
        for j in (i + 1)..=(per.len() - 1) {
            if per[i] > per[j] {
                count += 1;
            }
        }
    }
    if (count % 2) == 0 {
        result = false;
    }
    result
}
//формирует очеред-ное слагаемое в определителе
fn sum_matr_to_per(matrix: &Vec<Vec<f64>>, per: &Vec<usize>) -> f64 {
    let mut result = 1.0;
    for i in 0..per.len() {
        result *= matrix[i][per[i]];
    }
    if perestanovka(per) {
        result *= -1.0;
    }
    result
}
//рекурсивно формирует перестановки и ищет определитель/*
fn det_rec(matrix: &Vec<Vec<f64>>, per: &mut Vec<usize>, n: usize, n0: usize) -> f64 {
    let mut result = 0.0;
    for i in 0..n {
        if vkl(&per, i) {
            continue;
        }

        per[n0] = i;
        if (n - 1 == n0) {
            result = sum_matr_to_per(&matrix, &per);
        } else {
            result += det_rec(&matrix, per, n, n0 + 1);
        }
    }
    result
}
// подготавливает массив и запускает ре-курсию
//для нахождения определителя
fn det(matrix: &Vec<Vec<f64>>, n: usize) -> f64 {
    let mut per: Vec<usize> = vec![0; n];
    per[0] = 1;
    let result = det_rec(&matrix, &mut per, n, 0);
    result
}
//перемножение матриц
fn matrix_mult(
    first: &Vec<Vec<f64>>,
    second: &Vec<Vec<f64>>,
    n: usize,
    m: usize,
    k: usize,
) -> Vec<Vec<f64>> {
    let mut out_matrix = vec![vec![0.0; k]; n];
    for i in 0..n {
        for j in 0..k {
            let mut sum = 0.0;
            for l in 0..m {
                sum += first[i][l] * second[l][j];
            }
            out_matrix[i][j] = sum;
        }
    }
    out_matrix
}
// считает производную многочлена, переданного в массиве
fn proizvod(xar: &Vec<f64>) -> Vec<f64> {
    let mut proizv = vec![0.0; xar.len() - 1];
    for i in 0..(xar.len() - 2) {
        proizv[i] = xar[i] * ((xar.len() - i - 1) as f64);
    }
    proizv[xar.len() - 2] = xar[xar.len() - 2];
    proizv
}
// делит многочлен на одночлен (корень), тем самым уменьшая его степень
fn delenie(f: &Vec<f64>, koren: f64) -> Vec<f64> {
    let mut otv = vec![0.0; f.len()];
    otv[0] = f[0];
    for i in 1..f.len() {
        otv[i] = (koren * otv[i - 1]) + f[i];
    }
    otv
}
//подставляет число в многочлен
fn podstanovka(xar: &Vec<f64>, kor: f64) -> f64 {
    let mut result = 0.0;
    for i in 0..(xar.len() - 1) {
        result += kor.powi((xar.len() - 1 - i).try_into().unwrap()) * xar[i];
    }
    result += xar[xar.len() - 1];
    result
}
//находит решение многочлена
fn resh(mut xar: Vec<f64>) -> Vec<f64> {
    let mut dx: f64 = 0.0;
    let mut xn1 = 0.0;
    let mut otv = vec![0.0; xar.len() - 1];
    let mut f1 = vec![0.0; xar.len() - 1];
    for i in 1..xar.len() {
        let mut xn = 0.00001;
        f1 = proizvod(&xar);
        loop {
            dx = (-podstanovka(&xar, xn)) / podstanovka(&f1, xn);
            xn1 = dx + xn;
            xn = xn1;
            if dx.abs().le(&0.00001) {
                break;
            }
        }
        xar = delenie(&xar, xn1);
        let mut xar2 = vec![0.0; xar.len() - 1];
        for l in 0..xar2.len() {
            xar2[l] = xar[l];
        }
        xar = xar2.clone();
        otv[i - 1] = xn1;
    }
    otv
}
//находит значение очередного не-известного, считая сумму
//последующих элементов и деля её на элемент на главной диа-гонали
fn sum(matrix: &Vec<Vec<f64>>, mas: &Vec<f64>, p: usize) -> f64 {
    let mut result = 0.0;
    for i in (p + 1)..mas.len() {
        result += matrix[p][i] * mas[i];
    }
    result = (-result / matrix[p][p]);
    result
}
fn perest(matrix: &mut Vec<Vec<f64>>, p: usize, i: usize) -> bool {
    let mut tmp = 0.0;
    for u in (p + 1)..matrix.len() {
        if matrix[u][i] != 0.0 {
            for l in 0..matrix.len() {
                tmp = matrix[p][l];
                matrix[p][l] = matrix[u][l];
                matrix[u][l] = tmp;
            }
            return true;
        }
    }
    false
}
// делается проверка, если решение до этого было выбрано любое,
//а теперь выясняется что оно не подходит, то оно заменяет-ся 0
fn prov(matrix: &Vec<Vec<f64>>, b1: &mut Vec<f64>, k: usize, l: usize) {
    for i in (l + 1)..b1.len() {
        if matrix[k][i] != 0.0 {
            b1[i] = 0.0;
        }
    }
}
//приводим матрицу к ступенчатому виду и нахо-дим любое частное решение
fn stup(matrix: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let mut b = 0.0;
    for i in 0..(matrix.len() - 1) {
        for k in (i + 1)..matrix.len() {
            if matrix[i][i].abs() == 0.0 {
                if !perest(matrix, i, i) {
                    break;
                }
            }
            b = -matrix[k][i] / matrix[i][i];
            for j in 0..matrix.len() {
                matrix[k][j] = matrix[i][j] * b + matrix[k][j];
            }
            minim(matrix); //Возмодна проблема
        }
    }
    let mut b1 = vec![0.0; matrix.len()];
    for i in (0..=(matrix.len() - 1)).rev() {
        if matrix[i][i].abs() == 0.0 {
            b1[i] = 1.0;
            prov(matrix, &mut b1, i, i);
        } else {
            b1[i] = sum(&matrix, &b1, i);
        }
    }
    b1
}
fn minim(matrix: &mut Vec<Vec<f64>>) {
    for i in 0..matrix.len() {
        for j in 0..matrix.len() {
            if matrix[i][j].abs() < 0.0001 {
                matrix[i][j] = 0.0;
            }
        }
    }
}
//нормализуем массив
fn output(mut otv1: Vec<f64>) -> Vec<f64> {
    let mut s = 0.0;
    for i in 0..otv1.len() {
        s += otv1[i].powi(2);
    }
    s = s.sqrt();
    for i in 0..otv1.len() {
        otv1[i] /= s;
    }
    otv1.clone()
}
//находим собственные вектора для соб-ственных значений
fn sob_vect(matrix: &Vec<Vec<f64>>, otv: &Vec<f64>) -> Vec<Vec<f64>> {
    let mut out = Vec::new();
    let mut otv1 = vec![0.0; matrix.len() + 1];
    for k in 0..otv.len() {
        let mut matrix1 = matrix.clone();
        for i in 0..otv.len() {
            matrix1[i][i] -= otv[k];
        }
        minim(&mut matrix1);
        otv1 = stup(&mut matrix1);
        out.push(output(otv1))
    }
    out
}

fn laverre(matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut matrix1 = matrix.clone();
    let mut s = vec![0.0; matrix.len()];
    let mut p = vec![0.0; matrix.len() + 1];
    for i in 0..s.len() {
        if i != 0 {
            matrix1 = matrix_mult(&matrix1, &matrix, s.len(), s.len(), s.len());
        }
        for j in 0..s.len() {
            s[i] += matrix1[j][j];
        }
    }
    for i in 0..s.len() {
        p[i + 1] = s[i];
        if i != 0 {
            for j in 0..=(i - 1) {
                p[i + 1] = p[i + 1] + p[j + 1] * s[i - j - 1];
            }
        }
        p[i + 1] *= (-1.0 / ((i + 1) as f64));
    }
    p[0] = 1.0;
    p = resh(p);
    let sob_vectors = sob_vect(&matrix, &p);
    (p, sob_vectors)
}

fn fadeev(matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut matrix1 = matrix.clone();
    let mut matrix2 = matrix.clone();
    let mut s = vec![0.0; matrix[0].len()];
    let mut p = vec![0.0; matrix[0].len() + 1];
    for i in 0..matrix[0].len() {
        s[i] = 0.0;
        for j in 0..s.len() {
            s[i] = s[i] + matrix2[j][j];
        }
        s[i] = s[i] / (i as f64 + 1.0);
        matrix1 = matrix2.clone();
        for j in 0..matrix[0].len() {
            matrix1[j][j] = matrix2[j][j] - s[i];
        }
        matrix2 = matrix_mult(
            &matrix,
            &matrix1,
            matrix[0].len(),
            matrix[0].len(),
            matrix[0].len(),
        );
    }
    for i in 0..s.len() {
        p[i + 1] = -s[i];
    }
    p[0] = 1.0;
    p = resh(p);
    let sob_vectors = sob_vect(&matrix, &p);
    (p, sob_vectors)
}
// находим произведение матрицы на массив
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
//копируем и транспонируем матрицу, удаляя последнюю строчку
fn copy_trans(matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let mut result = vec![vec![0.0; matrix[0].len()]; matrix[0].len()];
    for i in 0..matrix[0].len() {
        for j in 0..matrix[0].len() {
            result[i][j] = matrix[j + 1][i];
        }
    }
    result
}
use rand::prelude::*;
fn krilov(matrix: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut p = vec![0.0; matrix.len() + 1];
    let mut l = 1.0;

    let mut a = vec![vec![0.0; matrix.len() - 1]; matrix.len() - 1];
    let mut y = vec![vec![0.0; matrix.len()]; matrix.len() + 1];
    let mut b = vec![0.0; matrix.len()];
    loop {
        a = matrix.clone();
        let mut rng = thread_rng();
        for j in 0..matrix.len() {
            loop {
                let rnd: i64 = rng.gen_range(0..100);
                y[matrix.len()][j] = (l * 10.0 * (rnd as f64 / 100.0)).round();
                if y[matrix.len()][j] >= 5.0 {
                    break;
                }
            }
        }
        l += 1.0;
        let mut tmp = vec![0.0; matrix.len()];
        for i in 1..=matrix.len() {
            for j in 0..matrix.len() {
                tmp[j] = y[matrix.len() - i + 1][j];
            }
            tmp = matr_umn_mas(&a, &tmp, matrix.len(), matrix.len(), matrix.len());
            for j in 0..matrix.len() {
                y[matrix.len() - i][j] = tmp[j];
            }
        }
        for j in 0..matrix.len() {
            b[j] = y[0][j];
        }
        a = copy_trans(&y);
        if det(&a, matrix.len()) == 0.0 {
            break;
        }
    }
    let mut a2 = vec![vec![0.0; matrix.len() + 1]; matrix.len()];
    for k in 0..matrix.len() {
        for j in 0..matrix.len() {
            a2[k][j] = a[k][j];
        }
    }
    a = a2.clone();
    b = methods::simple_gaus(a, b);
    for i in 1..p.len() {
        p[i] = -1.0 * b[i - 1];
    }
    p[0] = 1.0;
    p = resh(p);
    let sob_vectors = sob_vect(&matrix, &p);
    (p, sob_vectors)
}
