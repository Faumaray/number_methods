use std::io;
fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut free_elements = Vec::<f64>::new();
    let mut n = usize::MIN;
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
            n = v;
            }
            Err(..) => println!("Не подходящий формат: {}", trimmed),
            };
            for _i in 0..n
            {
                let mut temp = Vec::<f64>::new();
                for _j in 0..n
                {
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
    }
    else
    {
        matrix = vec![vec![3.72, 3.47,3.06], vec![4.47,4.10,3.63], vec![4.96,4.53,4.01]];
        free_elements = vec![30.74,36.80,40.79];
        n = 3;
    }
    println!("--------------------------------------------------------------------------------");
    println!("Исходные данные: ");
    println!("");
    for i in 0..n {
        for j in 0..n {
            if j != 0
            {
                if matrix[i][j] < 0.0
                {
                    print!(" - {}x{}", matrix[i][j].abs(),i);
                }
                else {
                    print!(" + {}x{}", matrix[i][j],i);
                }
            }
            else {
                print!(" {}x{}", matrix[i][j],i);
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
        let mut ans =0.0;
        for j in 0..n
        {
            ans += matrix[i][j]*answer_for_gauss[j];
        }
        println!("x[{}] = {}", i, ans-free_elements[i]);
    }
    println!("---------------------------------------------------------------------------------");
    println!("Корни системы по Холецкому: \n");
    for i in 0..n {
        println!("x[{}] = {}", i, answer_for_holeckiy[i]);
    }
    println!("Проверка");
    for i in 0..n {
        let mut ans =0.0;
        for j in 0..n
        {
            ans += matrix[i][j]*answer_for_holeckiy[j];
        }
        println!("x[{}] = {}", i, ans-free_elements[i]);
    }
}
fn holeckiy(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>) -> Vec<f64> {
    let n: usize = (&free_elements).len();
    let mut indent_matrix: Vec<Vec<f64>> = matrix.clone();
    let mut indent_free_elements: Vec<f64> = free_elements.clone();
    let mut c: Vec<Vec<f64>> = vec![vec![0.0_f64; n+1];n];
    let mut l: Vec<Vec<f64>> = vec![vec![0.0_f64; n];n];
    let mut y: Vec<f64> = vec![0.0_f64;n];
    let mut sum: f64 = 0.0_f64;
    //Умножение транспонированной матрицы на обычкую
    //Для создания симетричной положительно-определённой матрицы
    for i in 0..n {
        for j in 0..n {
            sum = 0.0_f64;
            for t in 0..n {
                sum = indent_matrix[t][j] * indent_matrix[t][i] + sum;
            }
            c[i][j] = sum;
        }
    }
    //{Тоже что и до этого только для правой части}
    for i in 0..n {
        for j in 0..n {
            y[i] = indent_matrix[j][i] * indent_free_elements[j] + y[i];
        }
    }
    
    for i in 0..n {
        for j in 0..n {
            indent_matrix[i][j] = c[i][j];
            indent_free_elements[i] = y[i];
        }
    }
    for i in 0..n {
        for j in 0..=i {
            sum = 0.0;
            if j != 0 {//Ветка для всех элементов матрицы l кроме 1 столбца
                for t in 0..=(j - 1) {
                    sum = sum + l[i][t] * l[j][t];
                }
                if i != j {
                    l[i][j] = (indent_matrix[i][j] - sum) / l[j][j];// для остальных кроме главной
                } else {
                    l[i][i] = (indent_matrix[i][i] - sum).sqrt();// для главной диагонали
                }
            } else {
                sum = sum + l[i][0] * l[j][0];
                if i != j {
                    l[i][j] = (indent_matrix[i][j] - sum) / l[j][j];// для остальных элементов первого столбца
                } else {
                    l[i][i] = (indent_matrix[i][i] - sum).sqrt();// для первого элемента
                }
            }
        }
    }
    for i in 0..n// Ly=b
    {
        sum = indent_free_elements[i];
        for j in 0..n
        {
            if i != j
            {
                sum -= l[i][j]*y[j];
            }
        }
        y[i] = sum/l[i][i];
    }
    for i in 0..n// Транспонирование l
    {
        for j in 0..n
        {
            l[i][j] = l[j][i];
        }
    }
    indent_free_elements = vec![0.0_f64;n]; // зануление свободных членов
    for i in (0..=(n-1)).rev() //Lt*x=y
    {
        sum = y[i];
        for j in 0..n
        {
            if i != j
            {
                sum -= l[i][j] * indent_free_elements[j];
            }
        }
        indent_free_elements[i] = sum/l[i][i];
    }
    return indent_free_elements.to_vec();
}
fn print(matrix: &Vec<Vec<f64>>)
{
    for row in matrix
    {
        println!("{:?}", row);
    }
}
fn gaus(mut orig: Vec<Vec<f64>>, free_elements: Vec<f64>) -> Vec<f64> {
    let n = (&free_elements).len();//Размерность начальной матрицы (строки)
    for i in 0..n 
    {
        orig[i].push(free_elements[i]);// Т.К. в решении этим методом свободные элементы учавствуют в преобразованиях
    }
    let mut matrix: Vec<Vec<f64>> = orig; 
    let mut matrix_clone = matrix.clone();//Матрица-дублер
    let mut x = vec![1.0;n];
    //Прямой ход (Зануление нижнего левого угла)
    for k in 0..n
    {
        for i in 0..=n
        {
            matrix_clone[k][i] /=matrix[k][k];//Деление k-строки на первый член !=0 для преобразования его в единицу
        }
        for i in (k+1)..n
        {
            let koef = matrix_clone[i][k]/matrix_clone[k][k];//Коэффициент умножения
            for j in 0..=n
            {
                matrix_clone[i][j]-= matrix_clone[k][j] * koef;//Зануление элементов матрицы ниже первого члена, преобразованного в единицу
            }
        }
        for i in 0..n
        {
            for j in 0..=n
            {
                matrix[i][j] = matrix_clone[i][j];
            }
        }
    }
    println!("После прямого хода:");
    print(&matrix);
     //Обратный ход (Зануление верхнего правого угла)
    for k in (0..=(n-1)).rev()
    {
        for i in (0..=n).rev()
        {
            matrix_clone[k][i] /= matrix[k][k];
        }
        if k == 0
        {
            continue;
        }
        for i in (0..=(k-1)).rev()
        {
            let koef = matrix_clone[i][k]/matrix_clone[k][k];
            for j in (0..=n).rev()
            {
                matrix_clone[i][j] -= matrix_clone[k][j] * koef;
            }
        }
    }
    println!("После обратного хода");
    print(&matrix_clone);
    for i in 0..n
    {
        x[i]=matrix_clone[i][n];
    }
    x
}
