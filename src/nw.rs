use std::cmp;

pub fn needleman_wunsch(str1: &str, str2: &str, gap_char: char) {
    let (best_candidate_sum_flag_matrix, alignment_score) = find_best_alignments(str1, str2);

    let aligned_strs = backtrace_best_alignment(str1, str2, gap_char, best_candidate_sum_flag_matrix);

    println!("\nbacktracing:");
    for string in aligned_strs {
        println!("{}", string);
    }

    println!("(alignment score = {})", alignment_score);
}

fn find_best_alignments(str1: &str, str2: &str) -> (Vec<Vec<u8>>, i32) {
    let mut score_path_matrix = vec![vec![0i32; str2.len()]; str1.len()];
    let mut best_candidate_sum_flag_matrix = vec![vec![0u8; str2.len()]; str1.len()];

    
    // parameter for score_path_matrix is to borrow per-call instead of for the entire scope where this closure exists
    let score_at = |row, col, score_path_matrix: &Vec<Vec<i32>>| match (row, col) {
        (-1i32, _) => -col - 1,
        (_, -1i32) => -row - 1,
        _ => score_path_matrix[row as usize][col as usize],
    };

    for (urow, char1) in str1.chars().enumerate() {
        for (ucol, char2) in str2.chars().enumerate() {
            let chars_match = char1 == char2;

            // Candidate sums
            
            let row = i32::try_from(urow).unwrap();
            let col = i32::try_from(ucol).unwrap();

            let match_mismatch_sum = score_at( row - 1, col - 1, &score_path_matrix)
                        + if chars_match {1} else {-1};
            let indel_row_sum = score_at(row - 1, col, &score_path_matrix) - 1;
            let indel_col_sum = score_at(row, col - 1, &score_path_matrix) - 1;

            
            let best_candidate_sum = cmp::max(match_mismatch_sum, cmp::max(indel_row_sum, indel_col_sum));

            score_path_matrix[urow][ucol] = best_candidate_sum;

            
            best_candidate_sum_flag_matrix[urow][ucol] = [match_mismatch_sum, indel_row_sum, indel_col_sum]
                    .iter()
                    .map(|&sum| sum == best_candidate_sum)
                    .enumerate()
                    .map(|(index, is_best)| (is_best as u8) << (index as u8))
                    .sum();
        }
    }

    println!("alignment score matrix:");
    for score_row in &score_path_matrix {
        println!("{:?}", score_row);
    }

    println!("\nbest candidate sum matrix:");
    for score_row in &best_candidate_sum_flag_matrix {
        println!("{:?}", score_row);
    }

    (best_candidate_sum_flag_matrix, *score_path_matrix.last().unwrap().last().unwrap())
}


fn backtrace_best_alignment(str1: &str, str2: &str, gap_char: char, best_candidate_sum_flag_matrix: Vec<Vec<u8>>) -> [String; 2] {
    // will be in reverse at first
    let mut aligned_str1 = Vec::new();
    let mut aligned_str2 = Vec::new();

    let mut row = str1.len() as i32 - 1;
    let mut col = str2.len() as i32 - 1;

    let mut n_remaining_chars1 = str1.len();
    let mut n_remaining_chars2 = str2.len();

    while row >= 0 && col >= 0 {
        let urow = usize::try_from(row).unwrap();
        let ucol = usize::try_from(col).unwrap();

        let best_candidate_sum_flag = best_candidate_sum_flag_matrix[urow][ucol];

        // Can pick any of the best candidate sum methods. Here we will just go with the first one found
        let best_candidate_sum_id = (0..3u8)
                .find(|&i| (best_candidate_sum_flag >> i) & 0b1 == 0b1)
                .unwrap();

        match best_candidate_sum_id {
            0 => {
                aligned_str1.push(str1.as_bytes()[urow] as char);
                aligned_str2.push(str2.as_bytes()[ucol] as char);

                row -= 1;
                col -= 1;
                n_remaining_chars1 -= 1;
                n_remaining_chars2 -= 1;
            },

            1 => {
                aligned_str1.push(str1.as_bytes()[urow] as char);
                aligned_str2.push(gap_char);

                row -= 1;
                n_remaining_chars1 -= 1;
            },

            2 => {
                aligned_str1.push(gap_char);
                aligned_str2.push(str2.as_bytes()[ucol] as char);

                col -= 1;
                n_remaining_chars2 -= 1;
            },

            _ => panic!(),
        };
    }

    let max_aligned_str_length = cmp::max(aligned_str1.len() + n_remaining_chars1, aligned_str2.len() + n_remaining_chars2);

    let aligned_str1: String =
            gap_char
            .to_string()
            .repeat(max_aligned_str_length - aligned_str1.len() - n_remaining_chars1)
            .chars()
            .chain((&str1[..n_remaining_chars1]).chars())
            .chain(aligned_str1.iter()
                    .rev()
                    .map(|&char| char)
            )
            .collect::<String>();

    let aligned_str2: String = 
            gap_char
            .to_string()
            .repeat(max_aligned_str_length - aligned_str2.len() - n_remaining_chars2)
            .chars()
            .chain((&str2[..n_remaining_chars2]).chars())
            .chain(aligned_str2.iter()
                    .rev()
                    .map(|&char| char)
            )
            .collect::<String>();

    [aligned_str1, aligned_str2]
}