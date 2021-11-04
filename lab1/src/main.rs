use std::io;
fn main() {
    let mut A = Vec::<Vec<f64>>::new();
    let mut B = Vec::<f64>::new();
    let mut N = usize::MIN;
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
    if inp == 0
    {
        println!("Введите размерность для квадратной матрицы");
        let mut input_text = String::new();
            io::stdin()
            .read_line(&mut input_text)
            .expect("Failed to read");
        let trimmed = input_text.trim();
        match trimmed.parse::<usize>() {
        Ok(v) => {
            N = v;
            }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
            };
            for _i in 0..N
            {
                let mut temp = Vec::<f64>::new();
                for _j in 0..N
                {
                    temp.push(0.0);
                }
                A.push(temp);
                B.push(0.0);
            }
        println!("Введите коэфы и свободные члены: ");
        for i in 0..N {
            for j in 0..N {
            let mut input_text = String::new();
            io::stdin()
            .read_line(&mut input_text)
            .expect("Failed to read");
        let trimmed = input_text.trim();
        match trimmed.parse::<f64>() {
        Ok(v) => {
            A[i][j] = v;
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
            B[i] = v;
        }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
        };
        }
    }
    else
    {
        A = vec![vec![3.72, 3.47,3.06], vec![4.47,4.10,3.63], vec![4.96,4.53,4.01]];
        B = vec![30.74,36.80,40.79];
        N = 3;
    }
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in 0..N {
        for j in 0..N {
            print!("A[{}][{}]={}  ||", i, j, A[i][j]);
        }
        print!("B[{}]={}", i, B[i]);
        println!("");
    }

    let answer_for_gauss = gaus(A.clone(), B.clone());
    let answer_for_holeckiy = holeckiy(A, B);
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Гауссу: \n");
    for i in 0..N {
        println!("x[{}] = {}", i, answer_for_gauss[i]);
    }
    println!("---------------------------------------------------------------------------------");
    println!("Корни системы по Холецкому: \n");
    for i in 0..N {
        println!("x[{}] = {}", i, answer_for_holeckiy[i]);
    }
}
fn holeckiy(A: Vec<Vec<f64>>, B: Vec<f64>) -> Vec<f64> {
    let N = (&B).len();
    let mut A1 = A.clone();
    let mut B1 = B.clone();
    let mut c = vec![vec![0.0; N+1];N];
    let mut L = vec![vec![0.0; N+1];N];
    let mut y = vec![0.0;N];
    let mut sum: f64 = 0.0;
    //Умножение матрицы на транспонированную
    //Для создания симетричной положительно-определённой матрицы
    for i in 0..N {
        for j in 0..N {
            sum = 0.0;
            for t in 0..N {
                sum = A1[t][j] * A1[t][i] + sum;
            }
            c[i][j] = sum;
        }
    }
    //{умножение правой части на транспонированную м-цу}
    for i in 0..N {
        for j in 0..N {
            y[i] = A1[j][i] * B1[j] + y[i];
        }
    }
    
    for i in 0..N {
        for j in 0..N {
            A1[i][j] = c[i][j];
            B1[i] = y[i];
        }
    }
    
    for i in 0..N {
        for j in 0..=i {
            sum = 0.0;
            if j != 0 {
                for t in 0..=(j - 1) {
                    sum = sum + L[i][t] * L[j][t];
                }
                if i != j {
                    L[i][j] = (A1[i][j] - sum) / L[j][j];
                } else {
                    L[i][i] = (A1[i][i] - sum).sqrt();
                }
            } else {
                sum = sum + L[i][0] * L[j][0];
                if i != j {
                    L[i][j] = (A1[i][j] - sum) / L[j][j];
                } else {
                    L[i][i] = (A1[i][i] - sum).sqrt();
                }
            }
        }
    }
    for i in 0..N {
        for j in 0..N {
            print!("L[{}][{}]={} ", i,j,L[i][j]);
        }
    }
    for i in 0..N {
        L[i][N] = B1[i];
    }
    B1[0] = L[0][N] / L[0][0];
    for i in 1..N {
        for j in 0..=(i - 1) {
            L[i][N] = L[i][N] - L[i][j] * B1[j];
            
        }
        B1[i] = L[i][N] / L[i][i];
        
    }
    for i in 0..N {
        for j in (i + 1)..N {
            L[i][j] = L[j][i];
            L[j][i] = 0.0;
        }
        L[i][N] = B1[i];
    }
    B1[N-1] = L[N-1][N] / L[N-1][N-1];
    
    for i in (0..=(N-1-1)).rev() {
        for j in (i + 1)..N {
            L[i][N] = L[i][N] - L[i][j] * B1[j];
        }
        B1[i] = L[i][N] / L[i][i];
    }
    return B1.to_vec();
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
