use methods::*;
use std::{io, vec};
fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut free_elements = Vec::<f64>::new();
    let mut size = usize::MIN;
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
                size = v;
            }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
        };
        for _i in 0..size {
            matrix.push(vec![0.0; size]);
            free_elements.push(0.0);
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
                    free_elements[i] = v;
                }
                Err(..) => println!("Не подходящий формат: {}", trimmed),
            };
        }
    } else {
        matrix = vec![vec![1.03, 0.993], vec![0.991, 0.051]];
        free_elements = vec![2.53, 2.43];
        size = 2;
    }
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in 0..size {
        for j in 0..size {
            if j != 0 {
                if matrix[i][j] < 0.0 {
                    print!(" - {}x{}", matrix[i][j].abs(), i);
                } else {
                    print!(" + {}x{}", matrix[i][j], i);
                }
            } else {
                print!(" {}x{}", matrix[i][j], i);
            }
        }
        print!(" = {}\n", free_elements[i]);
    }

    let answer_for_givens = givens(matrix.clone(), free_elements.clone(), size.clone());
    let answer_for_regular = regular(matrix.clone(), free_elements.clone(), size.clone());
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Гивенсу:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_givens[i]);
    }
    println!("Проверка");
    for i in 0..size {
        let mut ans = 0.0;
        for j in 0..size {
            ans += matrix[i][j] * answer_for_givens[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
    println!("\n---------------------------------------------------------------------------------");
    println!("Корни системы по Регулярному:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_regular[i]);
    }
    println!("Проверка");
    for i in 0..size {
        let mut ans = 0.0;
        for j in 0..size {
            ans += matrix[i][j] * answer_for_regular[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
}
