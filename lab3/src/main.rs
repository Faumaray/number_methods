use std::{io, vec};
fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut second = Vec::<f64>::new();
    let mut size = usize::MIN;
    println!("Ручной ввод?\n(1)Да;\n(2)Нет");
    let mut choice = String::new();
    io::stdin()
    .read_line(&mut choice)
    .expect("Failed to read");
    let mut inp = 0;
    match choice.trim().parse::<usize>() {
        Ok(v) => {
            inp = v;
        } ,
        Err(..) => eprintln!("Не подходящий формат: {}", choice.trim()),
    };
    if inp == 1
    {
        println!("Введите размерность для квадратной матрицы");
        let mut input_text = String::new();
            io::stdin()
            .read_line(&mut input_text)
            .expect("Failed to read");
        let trimmed = input_text.trim();
        match trimmed.parse::<usize>() {
        Ok(v) => {
            size = v;
            }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
            };
            for _i in 0..size
            {
                matrix.push(vec![0.0; size]);
                second.push(0.0);
            }
        println!("Введите коэфы и свободные члены: ");
        for i in 0..size {
            for j in 0..size {
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
            let mut input_text = String::new();
            io::stdin()
                .read_line(&mut input_text)
                .expect("Failed to read");
            let trimmed = input_text.trim();
            match trimmed.parse::<f64>() {
        Ok(v) => {
            second[i] = v;
        }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
        };
        }
    }
    else {
        matrix = vec![vec![1.03,0.993],vec![0.991,0.051]];
        second = vec![2.53,2.43];
        size = 2;
    }
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in 0..size {
        for j in 0..size {
            if j != 0
            {
                if matrix[i][j] < 0.0
                {
                    print!("{}*x{}", matrix[i][j],i);
                }
                else {
                    print!("+{}*x{}", matrix[i][j],i);
                }
            }
            else {
                print!("{}*x{}", matrix[i][j],i);
            }
        }
        print!("={}\n", second[i]);
    }

    let answer_for_givens = givens(matrix.clone(), second.clone(), size.clone());
    let answer_for_regular = regular(matrix.clone(), second.clone(), size.clone());
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Гивенсу:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_givens[i]);
    }
    println!("\n---------------------------------------------------------------------------------");
    println!("Корни системы по Регулярному:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_regular[i]);
    }
}
fn givens(matrix: Vec<Vec<f64>>, second: Vec<f64>, size: usize) -> Vec<f64>
{
    let mut answer = vec![0.0;size];
    let mut A = matrix.clone();
    let mut B = second.clone();
    let mut M = 0.0;
    let mut L;
    let mut R;
    for i in 0..(size-1)
    {
        for j in (i+1)..size
        {
            M = (matrix[i][i]*matrix[i][i]*matrix[j][i]*matrix[j][i]).sqrt();
            L = matrix[j][i] / M;// A12
            M = matrix[i][i]/ M;//B12
            for k in 0..size
            {
                R = matrix[i][k];
                A[i][k] = M * A[i][k] + L * A[j][k]; // a1j
                A[j][k] = M * A[j][k] - L * R;
            }
            R = B[i];
            B[i] = M * B[i] + L * B[j];
            B[j] = M * B[j] - L * R;
        }
    }
    println!("Матрица после вращения: {:?}", A);
    println!("B:{:?}", B);
    answer[size-1] = B[size-1] / A[size-1][size-1];
    for i in (0..=size-2).rev()
    {
        B[i] = second[i];
        A[i][i] = matrix[i][i];
        answer[i] = (B[i]-matrix[i][i+1] * answer[i+1])/ A[i][i];
    }
    answer
}
fn regular(matrix: Vec<Vec<f64>>, second: Vec<f64>, size: usize) -> Vec<f64>
{
    let mut answer = vec![0.0; size];
    const eps: f64 = 0.005;
    let mut a1 = vec![vec![0.0; size];size];
    let mut a2 = vec![vec![0.0; size];size];
    let mut b1 = vec![0.0;size];
    let mut b2 = vec![0.0;size];
    let mut x0 = vec![0.0;size];
    let mut s = 0.0;
    for i in  0..size
    {
        for k in 0..size
        {
            s = 0.0;
            for j in 0..size
            {
                s += matrix[j][i] * matrix[j][k];
            }
            a1[i][k] = s;
        }
    }
    for i in 0..size
    {
        s = 0.0;
        for j in 0..size
        {
            s += matrix[j][i] * second[j]
        }
        b1[i] = s;
    }
    let mut alfa = 0.0;
    b2 = vec![eps;size];
    let mut max = 0.0;
    loop {
        alfa += 0.00000001;
        a2 = a1.clone();
        for i in 0..size
        {
            a2[i][i] = a1[i][i] + alfa;
            b2[i] = b1[i] + alfa * x0[i];
        }
        a1 = a2.clone();
        b1 = b2.clone();
        b2 = gaus(a2.clone(), b2.clone());
        a2 = a1.clone();
        answer = b2.clone();
        x0 = answer.clone();
        b2 = b1.clone();
        b2 = gaus(a2, b2);
        max = (b2[1] - answer[1]).abs();
        for i in 1..size
        {
            if (b2[i] - answer[i]).abs() > max
            {
                max = (b2[i] - answer[i]).abs();
            }
        }
        if max <= eps
        {
            break;
        }
    }
    answer
}
fn gaus(mut A: Vec<Vec<f64>>, mut B: Vec<f64>) -> Vec<f64> {
    let N = (&B).len();
    let mut X = vec![1.0;N];
    //Вычитание из строк не другие строки умноженные на число
    for i in 0..N-1 {
        sort_rows(i, &mut A, &mut B, N);
        for j in (i + 1)..N {
            if A[i][i] != 0.0 {
                let mult_element = A[j][i] / A[i][i];
                for k in i..N {
                    A[j][k] = A[j][k] - (A[i][k] * mult_element);
                }
                B[j] = B[j] - (B[i] * mult_element);
            }
        }
    }
    for i in (0..=N-1).rev()
    {
        X[i] = B[i];
        for j in ((i+1)..N).rev()
        {
            X[i]=X[i]-(A[i][j]*X[j]);
        }
        X[i]=X[i]/A[i][i];
    }
    X
}
//Сортировка строки матрицы по возрастанию
fn sort_rows(sort_index: usize, matrix: &mut Vec<Vec<f64>>, right_part: &mut Vec<f64>, N:usize) {
    let mut max_element: f64 = matrix[sort_index][sort_index];
    let mut max_element_index = sort_index;
    for i in (sort_index + 1)..N {
        if matrix[i][sort_index] > max_element {
            max_element = matrix[i][sort_index];
            max_element_index = i;
        }
    }
    if max_element_index > sort_index {
        let mut temp: f64;
        temp = right_part[max_element_index];
        right_part[max_element_index] = right_part[sort_index];
        right_part[sort_index] = temp;
        for i in 0..N {
            temp = matrix[max_element_index][i];
            matrix[max_element_index][i] = matrix[sort_index][i];
            matrix[sort_index][i] = temp;
        }
    }
}
