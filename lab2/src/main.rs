use std::io;
use methods::*;
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
    if inp == 0 {
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
        matrix = vec![
            vec![16.0, 0.0688, 0.0829, 0.0970],
            vec![0.0496, 15.1, 0.0777, 0.0918],
            vec![0.0444, 0.0585, 14.2, 0.0867],
            vec![0.0393, 0.0534, 0.0674, 13.3],
        ];
        free_elements = vec![24.4781, 26.0849, 27.3281, 28.2078];
        size = 4;
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

    let answer_for_jacobi = jacobi(matrix.clone(), free_elements.clone());
    let answer_for_zejdel = zejdel(matrix.clone(), free_elements.clone());
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Якоби:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_jacobi.0[i]);
    }
    print!("Кол-во итераций = {}\n", answer_for_jacobi.1);
    println!("Проверка");
    for i in 0..size {
        let mut ans = 0.0;
        for j in 0..size {
            ans += matrix[i][j] * answer_for_jacobi.0[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
    println!("---------------------------------------------------------------------------------");
    println!("Корни системы по Зейделю:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_zejdel.0[i]);
    }
    print!("Кол-во итераций = {}\n", answer_for_zejdel.1);
    println!("Проверка");
    for i in 0..size {
        let mut ans = 0.0;
        for j in 0..size {
            ans += matrix[i][j] * answer_for_zejdel.0[j];
        }
        println!("x[{}] = {}", i, ans - free_elements[i]);
    }
    println!("---------------------------------------------------------------------------------");
    println!("Корни системы по обобщённому Зейделю:");
    for relax in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8].into_iter() {
        println!(
            "---------------------------------------------------------------------------------"
        );
        println!("Ответ при релаксации = {}", &relax);
        let answer_for_zejdelex = zejdelex(matrix.clone(), free_elements.clone(), relax.clone());

        for i in 0..size {
            print!("x{} = {} ", i, answer_for_zejdelex.0[i]);
        }
        print!("Кол-во итераций = {} w={}\n", answer_for_zejdelex.1, relax);
        println!("Проверка");
        for i in 0..size {
            let mut ans = 0.0;
            for j in 0..size {
                ans += matrix[i][j] * answer_for_zejdelex.0[j];
            }
            println!("x[{}] = {}", i, ans - free_elements[i]);
        }
    }
}
