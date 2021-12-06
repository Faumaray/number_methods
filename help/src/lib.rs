pub static EPS: f64 = 0.0001;
// LAB1
pub fn holeckiy(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>) -> Vec<f64> {
    let n: usize = (&free_elements).len();
    let mut indent_matrix: Vec<Vec<f64>> = matrix.clone();
    let mut indent_free_elements: Vec<f64> = free_elements.clone();
    let mut c: Vec<Vec<f64>> = vec![vec![0.0_f64; n + 1]; n];
    let mut l: Vec<Vec<f64>> = vec![vec![0.0_f64; n]; n];
    let mut y: Vec<f64> = vec![0.0_f64; n];
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
            if j != 0 {
                //Ветка для всех элементов матрицы l кроме 1 столбца
                for t in 0..=(j - 1) {
                    sum = sum + l[i][t] * l[j][t];
                }
                if i != j {
                    l[i][j] = (indent_matrix[i][j] - sum) / l[j][j]; // для остальных кроме главной
                } else {
                    l[i][i] = (indent_matrix[i][i] - sum).sqrt(); // для главной диагонали
                }
            } else {
                sum = sum + l[i][0] * l[j][0];
                if i != j {
                    l[i][j] = (indent_matrix[i][j] - sum) / l[j][j]; // для остальных элементов первого столбца
                } else {
                    l[i][i] = (indent_matrix[i][i] - sum).sqrt(); // для первого элемента
                }
            }
        }
    }
    for i in 0..n
    // Ly=b
    {
        sum = indent_free_elements[i];
        for j in 0..n {
            if i != j {
                sum -= l[i][j] * y[j];
            }
        }
        y[i] = sum / l[i][i];
    }
    for i in 0..n
    // Транспонирование l
    {
        for j in 0..n {
            l[i][j] = l[j][i];
        }
    }
    indent_free_elements = vec![0.0_f64; n]; // зануление свободных членов
    for i in (0..=(n - 1)).rev()
    //Lt*x=y
    {
        sum = y[i];
        for j in 0..n {
            if i != j {
                sum -= l[i][j] * indent_free_elements[j];
            }
        }
        indent_free_elements[i] = sum / l[i][i];
    }
    return indent_free_elements.to_vec();
}
fn print(matrix: &Vec<Vec<f64>>) {
    for row in matrix {
        println!("{:?}", row);
    }
}
pub fn gaus(mut orig: Vec<Vec<f64>>, free_elements: Vec<f64>) -> Vec<f64> {
    let n = (&free_elements).len(); //Размерность начальной матрицы (строки)
    for i in 0..n {
        orig[i].push(free_elements[i]); // Т.К. в решении этим методом свободные элементы учавствуют в преобразованиях
    }
    let mut matrix: Vec<Vec<f64>> = orig;
    let mut matrix_clone = matrix.clone(); //Матрица-дублер
    let mut x = vec![1.0; n];
    //Прямой ход (Зануление нижнего левого угла)
    for k in 0..n {
        for i in 0..=n {
            matrix_clone[k][i] /= matrix[k][k]; //Деление k-строки на первый член !=0 для преобразования его в единицу
        }
        for i in (k + 1)..n {
            let koef = matrix_clone[i][k] / matrix_clone[k][k]; //Коэффициент умножения
            for j in 0..=n {
                matrix_clone[i][j] -= matrix_clone[k][j] * koef; //Зануление элементов матрицы ниже первого члена, преобразованного в единицу
            }
        }
        for i in 0..n {
            for j in 0..=n {
                matrix[i][j] = matrix_clone[i][j];
            }
        }
    }
    println!("После прямого хода:");
    print(&matrix);
    //Обратный ход (Зануление верхнего правого угла)
    for k in (0..=(n - 1)).rev() {
        for i in (0..=n).rev() {
            matrix_clone[k][i] /= matrix[k][k];
        }
        if k == 0 {
            continue;
        }
        for i in (0..=(k - 1)).rev() {
            let koef = matrix_clone[i][k] / matrix_clone[k][k];
            for j in (0..=n).rev() {
                matrix_clone[i][j] -= matrix_clone[k][j] * koef;
            }
        }
    }
    println!("После обратного хода");
    print(&matrix_clone);
    for i in 0..n {
        x[i] = matrix_clone[i][n];
    }
    x
}

// LAB2
pub fn zejdelex(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>, relax: f64) -> (Vec<f64>, usize) // пересмотреть решает вроде верно, но не похоже на правильный алгоритм
{
    let size: usize = free_elements.len();
    let mut answer = vec![0.0; size];
    let mut prev = vec![0.0; size];
    let mut iter_count = 0;
    let mut norm = f64::MAX;
    while norm >= EPS {
        iter_count += 1;
        for i in 0..size {
            answer[i] = free_elements[i];
            for j in 0..size {
                if i != j {
                    answer[i] = answer[i] - matrix[i][j] * answer[j];
                }
            }
            answer[i] /= matrix[i][i];
            answer[i] = relax * answer[i] + (1.0 - relax) * prev[i];
        }
        for i in 0..size {
            if (answer[i] - prev[i]).abs() < norm {
                norm = (answer[i] - prev[i]).abs();
            }
            prev[i] = answer[i];
        }
    }
    (answer, iter_count)
}
pub fn zejdel(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>) -> (Vec<f64>, usize) {
    let size: usize = free_elements.len();
    let mut answer = vec![0.0; size];
    let mut prev = vec![0.0; size];
    let mut iter_count = 0;
    let mut norm = f64::MAX;
    while norm.sqrt() >= EPS {
        norm = 0.0;
        iter_count += 1;
        for i in 0..size {
            prev[i] = answer[i];
        }
        for i in 0..size {
            let mut var = 0.0;
            for j in 0..size {
                if j != i {
                    var += matrix[i][j] * answer[j];
                }
            }
            answer[i] = (free_elements[i] - var) / matrix[i][i];
            norm += (answer[i] - prev[i]).powi(2);
        }
    }
    (answer, iter_count)
}
pub fn jacobi(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>) -> (Vec<f64>, usize) {
    let size: usize = free_elements.len();
    let mut answer: Vec<f64> = vec![1.0; size];
    let mut temp_answer: Vec<f64> = vec![0.0; size];
    let mut norm = f64::MAX;
    let mut iter_count = 0;
    while norm >= EPS {
        iter_count += 1;
        for i in 0..size {
            temp_answer[i] = 0.0 - free_elements[i];
            for g in 0..size {
                if i != g {
                    temp_answer[i] += matrix[i][g] * answer[g];
                }
            }
            temp_answer[i] /= 0.0 - matrix[i][i];
        }
        norm = (answer[0] - temp_answer[0]).abs();
        for h in 0..size {
            if (answer[h] - temp_answer[h]).abs() > norm {
                norm = (answer[h] - temp_answer[h]).abs();
            }
            answer[h] = temp_answer[h];
        }
    }
    (answer, iter_count)
}
// LAB3
pub fn givens(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>, size: usize) -> Vec<f64> {
    let mut answer = vec![0.0; size];
    let mut indent_matrix = matrix.clone();
    let mut indent_free_elements = free_elements.clone();
    let mut l;
    let mut r;
    for i in 0..(size - 1) {
        for j in (i + 1)..size {
            let mut m = (matrix[i][i] * matrix[i][i] * matrix[j][i] * matrix[j][i]).sqrt();
            l = matrix[j][i] / m; // A12
            m = matrix[i][i] / m; //B12
            for k in 0..size {
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
            if j != 0 {
                if indent_matrix[i][j] < 0.0 {
                    print!(" - {}", indent_matrix[i][j].abs());
                } else {
                    print!(" + {}", indent_matrix[i][j]);
                }
            } else {
                print!("{}", indent_matrix[i][j]);
            }
        }
        print!(" = {}\n", indent_free_elements[i]);
    }
    answer[size - 1] = indent_free_elements[size - 1] / indent_matrix[size - 1][size - 1];
    for i in (0..=size - 2).rev() {
        indent_free_elements[i] = free_elements[i];
        indent_matrix[i][i] = matrix[i][i];
        answer[i] =
            (indent_free_elements[i] - matrix[i][i + 1] * answer[i + 1]) / indent_matrix[i][i];
    }
    answer
}
pub fn regular(matrix: Vec<Vec<f64>>, free_elements: Vec<f64>, size: usize) -> Vec<f64> {
    let mut answer = vec![0.0; size];
    let mut a1 = vec![vec![0.0; size]; size];
    let mut b1 = vec![0.0; size];
    let mut x0 = vec![0.0; size];
    for i in 0..size {
        for k in 0..size {
            let mut s = 0.0;
            for j in 0..size {
                s += matrix[j][i] * matrix[j][k];
            }
            a1[i][k] = s;
        }
    }
    for i in 0..size {
        let mut s = 0.0;
        for j in 0..size {
            s += matrix[j][i] * free_elements[j]
        }
        b1[i] = s;
    }
    let mut alfa = 0.0;
    let mut b2 = vec![EPS; size];
    let mut max = f64::MAX;
    while max >= EPS {
        alfa += 0.00000001;
        let mut a2 = a1.clone();
        for i in 0..size {
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
        for i in 1..size {
            if (b2[i] - answer[i]).abs() > max {
                max = (b2[i] - answer[i]).abs();
            }
        }
    }
    answer
}

// Alternative
pub fn simple_gaus(
    mut indent_matrix: Vec<Vec<f64>>,
    mut indent_free_elements: Vec<f64>,
) -> Vec<f64> {
    let size = (&indent_free_elements).len();
    let mut answer = vec![1.0; size];
    //Вычитание из строк не другие строки умноженные на число
    for i in 0..size - 1 {
        sort_rows(i, &mut indent_matrix, &mut indent_free_elements, size);
        for j in (i + 1)..size {
            if indent_matrix[i][i] != 0.0 {
                let mult_element = indent_matrix[j][i] / indent_matrix[i][i];
                for k in i..size {
                    indent_matrix[j][k] =
                        indent_matrix[j][k] - (indent_matrix[i][k] * mult_element);
                }
                indent_free_elements[j] =
                    indent_free_elements[j] - (indent_free_elements[i] * mult_element);
            }
        }
    }
    for i in (0..=size - 1).rev() {
        answer[i] = indent_free_elements[i];
        for j in ((i + 1)..size).rev() {
            answer[i] = answer[i] - (indent_matrix[i][j] * answer[j]);
        }
        answer[i] = answer[i] / indent_matrix[i][i];
    }
    answer
}
//Сортировка строки матрицы по возрастанию
fn sort_rows(
    sort_index: usize,
    matrix: &mut Vec<Vec<f64>>,
    right_part: &mut Vec<f64>,
    size: usize,
) {
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
pub fn f1(y1: f64, y2: f64, x: f64) -> f64 {
    (y1 + x) / (y1.powi(2) + y2.powi(2))
}
pub fn f2(y1: f64, y2: f64, x: f64) -> f64 {
    (y1 + x * y2).cos()
}

// Lab 8
pub fn eiler(
    a: f64,
    b: f64,
    n: usize,
    kolfun: usize,
    mut x: f64,
    f: Vec<&dyn Fn(f64, f64, f64) -> f64>,
    y_1: &mut Vec<Vec<f64>>,
) -> Vec<(f64, f64, f64)> {
    let mut out: Vec<(f64, f64, f64)> = Vec::new();
    let t: f64 = (b - a) / (n as f64);
    for i in 1..n + 1 {
        for k in 0..kolfun {
            let value = y_1[k][i - 1] + t * f[k](x, y_1[0][i - 1], y_1[1][i - 1]);
            y_1[k].push(value);
        }
        x += t;
        out.push((x, y_1[0][i], y_1[1][i]));
    }
    out.push((a, y_1[0][1], y_1[1][2]));
    out
}
pub fn prognoz(
    a: f64,
    b: f64,
    n: usize,
    kolfun: usize,
    mut x: f64,
    f: Vec<&dyn Fn(f64, f64, f64) -> f64>,
    y_1: &mut Vec<Vec<f64>>,
) -> Vec<(f64, f64, f64)> {
    let mut out: Vec<(f64, f64, f64)> = Vec::new();
    let t: f64 = (b - a) / (n as f64);
    for i in 1..n + 1 {
        for k in 0..kolfun {
            let value = y_1[k][i - 1]
                + t * f[k](
                    x + 0.5 * t,
                    y_1[0][i - 1] + 0.5 * t * y_1[0][i - 1],
                    y_1[1][i - 1] + 0.5 * t * y_1[1][i - 1],
                );
            y_1[k].push(value);
        }
        x += t;
        out.push((x, y_1[0][i], y_1[1][i]));
    }
    out.push((a, y_1[0][1], y_1[1][2]));
    out
}
pub fn runge_kut(
    a: f64,
    h: f64,
    n: usize,
    kolfun: usize,
    mut x: f64,
    f: Vec<&dyn Fn(f64, f64, f64) -> f64>,
    y_1: &mut Vec<Vec<f64>>,
) -> Vec<(f64, f64, f64)> {
    let mut k = [0.0; 4];
    let mut out: Vec<(f64, f64, f64)> = Vec::new();
    for i in 1..n + 1 {
        for j in 0..kolfun {
            k[0] = f[j](x, y_1[0][i - 1], y_1[1][i - 1]);
            k[1] = f[j](
                x + h / 2.0,
                y_1[0][i - 1] + h / 2.0 * k[1],
                y_1[1][i - 1] + h / 2.0 * k[1],
            );
            k[2] = f[j](
                x + h / 2.0,
                y_1[0][i - 1] + h / 2.0 * k[2],
                y_1[1][i - 1] + h / 2.0 * k[2],
            );
            k[3] = f[j](x + h, y_1[0][i - 1] + h * k[3], y_1[1][i - 1] + h * k[3]);
            let value = y_1[j][i - 1] + h / 6.0 * (k[0] + 2.0 * k[1] + 2.0 * k[2] + k[3]);
            y_1[j].push(value);
        }
        x += h;
        out.push((x, y_1[0][i], y_1[1][i]));
    }
    out.push((a, y_1[0][1], y_1[1][2]));
    out
}
