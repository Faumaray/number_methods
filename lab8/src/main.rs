use methods::{eiler, f1, f2, prognoz, runge_kut};

fn main() {
    let func_list: Vec<&dyn Fn(f64, f64, f64) -> f64> = vec![&f1, &f2];
    let answer_eiler = eiler(
        2.0,
        4.0,
        10,
        2,
        2.0,
        func_list.clone(),
        &mut vec![vec![1.0], vec![-1.0]],
    );
    let answer_prognoz = prognoz(
        2.0,
        4.0,
        10,
        2,
        2.0,
        func_list.clone(),
        &mut vec![vec![1.0], vec![-1.0]],
    );
    let answer_runge_kut = runge_kut(
        2.0,
        4.0,
        10,
        2,
        2.0,
        func_list.clone(),
        &mut vec![vec![1.0], vec![-1.0]],
    );
    print_answer(answer_eiler, "Эйлеру");
    print_answer(answer_prognoz, "Прогнозу");
    print_answer(answer_runge_kut, "Рунге");
}
fn print_answer(answer: Vec<(f64, f64, f64)>, name: &str) {
    println!("-----------------------------------------------------------------");
    println!("Ответы по {}", name);
    for i in 0..answer.len() {
        println!(
            "x: {:.02}, y1: {:.05}, y2: {:.05}",
            answer[i].0, answer[i].1, answer[i].2
        );
    }
}
