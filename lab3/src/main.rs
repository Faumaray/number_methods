use std::{io, vec};
fn main() {
    let mut matrix = Vec::<Vec<f64>>::new();
    let mut free_elements = Vec::<f64>::new();
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
    }
    else {
        matrix = vec![vec![1.03,0.993],vec![0.991,0.051]];
        free_elements = vec![2.53,2.43];
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

    let answer_for_givens = givens(matrix.clone(), free_elements.clone(), size.clone());
    let answer_for_regular = regular(matrix.clone(), free_elements.clone(), size.clone());
    println!("--------------------------------------------------------------------------------");
    println!("Корни системы по Гивенсу:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_givens[i]);
    }
    println!("Проверка");
        for i in 0..size {
            let mut ans =0.0;
            for j in 0..size
            {
                ans += matrix[i][j]*answer_for_givens[j];
            }
            println!("x[{}] = {}", i, ans-free_elements[i]);
        }
    println!("\n---------------------------------------------------------------------------------");
    println!("Корни системы по Регулярному:");
    for i in 0..size {
        print!("x{} = {} ", i, answer_for_regular[i]);
    }
    println!("Проверка");
        for i in 0..size {
            let mut ans =0.0;
            for j in 0..size
            {
                ans += matrix[i][j]*answer_for_regular[j];
            }
            println!("x[{}] = {}", i, ans-free_elements[i]);
        }
}
fn givens(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>, size: usize) -> Vec<f64>
{
    let mut answer = vec![0.0;size];
    let mut indent_matrix = matrix.clone();
    let mut indent_free_elements = free_elements.clone();
    let mut l;
    let mut r;
    for i in 0..(size-1)
    {
        for j in (i+1)..size
        {
            let mut m = (matrix[i][i]*matrix[i][i]*matrix[j][i]*matrix[j][i]).sqrt();
            l = matrix[j][i] / m;// A12
            m = matrix[i][i]/ m;//B12
            for k in 0..size
            {
                r = matrix[i][k];
                indent_matrix[i][k] = m * indent_matrix[i][k] + l * indent_matrix[j][k]; // a1j
                indent_matrix[j][k] = m * indent_matrix[j][k] - l * r;
            }
            r = indent_free_elements[i];
            indent_free_elements[i] = m * indent_free_elements[i] + l * indent_free_elements[j];
            indent_free_elements[j] = m * indent_free_elements[j] - l * r;
        }
    }
    println!("Данные после вращения: \n");
    for i in 0..size {
        for j in 0..size {
            if j != 0
            {
                if indent_matrix[i][j] < 0.0
                {
                    print!(" - {}", indent_matrix[i][j].abs());
                }
                else {
                    print!(" + {}", indent_matrix[i][j]);
                }
            }
            else {
                print!("{}", indent_matrix[i][j]);
            }
        }
        print!(" = {}\n", indent_free_elements[i]);
    }
    answer[size-1] = indent_free_elements[size-1] / indent_matrix[size-1][size-1];
    for i in (0..=size-2).rev()
    {
        indent_free_elements[i] = free_elements[i];
        indent_matrix[i][i] = matrix[i][i];
        answer[i] = (indent_free_elements[i]-matrix[i][i+1] * answer[i+1])/ indent_matrix[i][i];
    }
    answer
}
fn regular(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>, size: usize) -> Vec<f64>
{
    let mut answer = vec![0.0; size];
    const CUSTOM_EPS: f64 = 0.005;
    let mut a1 = vec![vec![0.0; size];size];
    let mut b1 = vec![0.0;size];
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
            s += matrix[j][i] * free_elements[j]
        }
        b1[i] = s;
    }
    let mut alfa = 0.0;
    let mut b2 = vec![CUSTOM_EPS;size];
    let mut max = f64::MAX;
    while max>=CUSTOM_EPS {
        alfa += 0.00000001;
        let mut a2 = a1.clone();
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

    }
    answer
}
fn gaus(mut indent_matrix: Vec<Vec<f64>>, mut indent_free_elements: Vec<f64>) -> Vec<f64> {
    let size = (&indent_free_elements).len();
    let mut answer = vec![1.0;size];
    //Вычитание из строк не другие строки умноженные на число
    for i in 0..size-1 {
        sort_rows(i, &mut indent_matrix, &mut indent_free_elements, size);
        for j in (i + 1)..size {
            if indent_matrix[i][i] != 0.0 {
                let mult_element = indent_matrix[j][i] / indent_matrix[i][i];
                for k in i..size {
                    indent_matrix[j][k] = indent_matrix[j][k] - (indent_matrix[i][k] * mult_element);
                }
                indent_free_elements[j] = indent_free_elements[j] - (indent_free_elements[i] * mult_element);
            }
        }
    }
    for i in (0..=size-1).rev()
    {
        answer[i] = indent_free_elements[i];
        for j in ((i+1)..size).rev()
        {
            answer[i]=answer[i]-(indent_matrix[i][j]*answer[j]);
        }
        answer[i]=answer[i]/indent_matrix[i][i];
    }
    answer
}
//Сортировка строки матрицы по возрастанию
fn sort_rows(sort_index: usize, matrix: &mut Vec<Vec<f64>>, right_part: &mut Vec<f64>, size:usize) {
    let mut max_element: f64 = matrix[sort_index][sort_index];
    let mut max_element_index = sort_index;
    for i in (sort_index + 1)..size {
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
        for i in 0..size {
            temp = matrix[max_element_index][i];
            matrix[max_element_index][i] = matrix[sort_index][i];
            matrix[sort_index][i] = temp;
        }
    }
}
