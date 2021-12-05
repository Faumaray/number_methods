use methods::*;
use std::io;
fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut free_elements = Vec::<f64>::new();
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
            free_elements.push(0.0);
        }
        println!("Введите коэфы и свободные члены: ");
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
        matrix = vec![
            vec![3.72, 3.47, 3.06],
            vec![4.47, 4.10, 3.63],
            vec![4.96, 4.53, 4.01],
        ];
        free_elements = vec![30.74, 36.80, 40.79];
        n = 3;
    }
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in 0..n {
        for j in 0..n {
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

    let answer_for_gauss = gaus(matrix.clone(), free_elements.clone());
    let answer_for_holeckiy = holeckiy(matrix.clone(), free_elements.clone());
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Гауссу: \n");
    for i in 0..n {
        println!("x[{}] = {}", i, answer_for_gauss[i]);
    }
    println!("Проверка");
    for i in 0..n {
        let mut ans = 0.0;
        for j in 0..n {
            ans += matrix[i][j] * answer_for_gauss[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
    println!("---------------------------------------------------------------------------------");
    println!("Корни системы по Холецкому: \n");
    for i in 0..n {
        println!("x[{}] = {}", i, answer_for_holeckiy[i]);
    }
    println!("Проверка");
    for i in 0..n {
        let mut ans = 0.0;
        for j in 0..n {
            ans += matrix[i][j] * answer_for_holeckiy[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
}
