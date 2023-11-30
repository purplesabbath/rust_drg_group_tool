// use serde::Deserialize;
// use serde_json::Value;
use polars::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::io::{self, Write};
use std::path::Path;
use std::string::String;
// use serde::de::Unexpected::Str;
// use serde_json::Value::String;

// 读取分组方案=======================================================================================
// 读取JSON文件为字典(HashMap), 键为MDC编码, 值为MDC下的主诊断HashSet
fn read_file_as_str_to_set<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, HashSet<String>>, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    Ok(u)
}

// 读取JSON文件为字典(HashMap), 键为MDC编码, 值为string
fn read_file_as_str_to_str<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    Ok(u)
}

// 读取JSON文件为字典(HashMap), 键为MDC编码, 值为向量
fn read_file_as_str_to_tuple<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, Vec<String>>, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `User`.
    let u = serde_json::from_reader(reader)?;

    Ok(u)
}

// 读取所有手术列表
fn read_icd9_to_vec<P: AsRef<Path>>(file_path: P) -> Result<HashSet<String>, Box<dyn Error>> {
    let contents = fs::read_to_string(file_path)?;
    let vec: HashSet<String> = contents.split(',').map(|s| s.to_string()).collect();
    Ok(vec)
}

// 病例结构===========================================================================================
struct DrgCase {
    id: String,               // 病例ID
    main_dis: String,         // 主诊断编码(必填)
    main_opt: String,         // 主手术编码(手术病例必填)
    other_dis: Vec<String>,   // 其他诊断编码(列表)
    other_opt: Vec<String>,   // 其他手术编码(列表)
    sex: i32,                 // 性别(0 => 女, 1 => 男)
    age: f64,                 // 年龄(不足一岁以小于1小数表示, 出生天数/365)
    weight: i32,              // 体重
    all_dis: HashSet<String>, // 所有的诊断
    all_opt: HashSet<String>, // 所有的手术
}

impl DrgCase {
    // 初始化方法
    fn new(
        admission_number: String,
        principal_diagnosis: String,
        principal_operation: String,
        other_diagnosis: Vec<String>,
        other_operation: Vec<String>,
        gender: i32,
        old: f64,
        mass: i32,
    ) -> Self {
        Self {
            id: admission_number,
            main_dis: principal_diagnosis,
            main_opt: principal_operation,
            other_dis: other_diagnosis,
            other_opt: other_operation,
            sex: gender,
            age: old,
            weight: mass,
            all_dis: HashSet::new(), // 初始化为空
            all_opt: HashSet::new(), // 初始化为空
        }
    }

    // 检查病例是否有主手术
    fn no_surgery(&self) -> bool {
        return self.main_opt == "";
    }

    // 检查病例是否有其他手术
    fn no_other_surgery(&self) -> bool {
        return self.other_opt.len() == 0;
    }

    // 检查病例是否有其他诊断
    fn no_other_diagnosis(&self) -> bool {
        return self.other_dis.len() == 0;
    }

    // 合并主诊断与其他诊断为一个set
    fn concat_dis(&mut self) {
        let mut temp_dis = self.other_dis.clone();
        let principle_dis = self.main_dis.clone();
        self.all_dis.insert(principle_dis);
        temp_dis.retain(|r| self.all_dis.insert(r.to_string()))
    }

    // 合并主诊断与其他诊断为一个set
    fn concat_opt(&mut self) {
        let mut temp_opt = self.other_opt.clone();
        match self.main_opt.as_str() {
            "" => (),
            _ => {
                let mut temp_opt = self.other_opt.clone();
                let principle_opt = self.main_opt.clone();
                self.all_opt.insert(principle_opt);
                temp_opt.retain(|r| self.all_opt.insert(r.to_string()))
            }
        };
    }
}

// 分组逻辑===========================================================================================
// MDC判断>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// MDCA先期分组
fn is_mdca(
    record: &DrgCase,                                // 病例
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    all_opt_list: &HashSet<String>,                  // 全部手术列表
    adrg_type_dict: &HashMap<String, Vec<String>>,   // ADRG类型及对应入组类型
    mdc_name: String,
) -> String {
    // 如果病例没有手术不符合入组条件
    if record.no_surgery() {
        return String::from("KBBZ");
    }
    let mut pred = String::from("KBBZ");
    let adrg_list = vec![
        "AA1", "AB1", "AC1", "AD1", "AE1", "AF1", "AG1", "AG2", "AH1",
    ];
    for adrg in adrg_list {
        pred = process_adrg(
            record,
            adrg_dis_opt,
            all_opt_list,
            adrg_type_dict,
            adrg.to_string(),
        );
        if pred != "KBBZ" {
            break;
        }
    }
    return pred;
}

// MDCZ多发创伤
fn is_mdcz(
    record: &DrgCase,
    mdcz_adrg_dis_dict: &HashMap<String, HashSet<String>>,
    mdc_name: String,
) -> String {
    // 如果病例没有其他诊断不符合入组条件(需要至少两个诊断)
    if record.no_other_diagnosis() {
        return String::from("KBBZ");
    }

    // 需要主诊断或其他诊断分别在两个不同部位
    let mut condition_conut = 0;
    for frag in [
        "head_dis",
        "chest_dis",
        "belly_dis",
        "urinary_dis",
        "reproductive_dis",
        "torso_spine_dis",
        "upper_limb_dis",
        "lower_limb_dis",
        "bon_dis",
    ] {
        if !(mdcz_adrg_dis_dict[frag].is_disjoint(&record.all_dis)) {
            condition_conut += 1;
        }
    }
    return if condition_conut >= 1 {
        String::from("MDCZ")
    } else {
        String::from("KBBZ")
    };
}

// MDCP需要根据的年龄进行判断的
fn is_age_mdc(
    record: &DrgCase,
    mdc_dis: &HashMap<String, HashSet<String>>,
    mdc_name: String,
) -> String {
    if (record.age <= (29 / 365) as f64) & (mdc_dis[&mdc_name].contains(&record.main_dis)) {
        return mdc_name;
    }
    return String::from("KBBZ");
}

// MDCM和MDCN需要根据性别进行判断的
fn is_sex_mdc(
    record: &DrgCase,
    mdc_dis: &HashMap<String, HashSet<String>>,
    mdc_name: String,
) -> String {
    if (mdc_name == "MDCM") & (record.sex == 1) & (mdc_dis[&mdc_name].contains(&record.main_dis)) {
        return mdc_name;
    } else if (mdc_name == "MDCN")
        & (record.sex == 0)
        & (mdc_dis[&mdc_name].contains(&record.main_dis))
    {
        return mdc_name;
    } else {
        return String::from("KBBZ");
    }
}

// 普通的根据主诊断入组的MDC
fn is_common_mdc(
    record: &DrgCase,
    mdc_dis: &HashMap<String, HashSet<String>>,
    mdc_name: String,
) -> String {
    return if mdc_dis[&mdc_name].contains(&record.main_dis) {
        mdc_name
    } else {
        String::from("KBBZ")
    };
}

// ADRG判断>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// 只根据主手术入组的ADRG
fn is_common_surgery_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 如果无手术, 则无法入组
    if record.no_surgery() {
        return "KBBZ".to_string();
    }
    // 主手术在该ADRG的主手术表中
    let verb: String = adrg_name.to_string() + "_opt"; // 字符串拼接生成键
    return if adrg_dis_opt_dict[&verb].contains(&record.main_opt) {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 只根据主诊断入组的ADRG
fn is_common_diagnosis_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 主诊断在该ADRG的主诊断表中
    let verb = adrg_name.to_string() + "_dis";
    return if adrg_dis_opt_dict[&verb].contains(&record.main_dis) {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组的ADRG => 主诊断 + 主手术
fn is_both_mdis_and_mopt_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 如果病例没有主手术, 无需判断
    if record.no_surgery() {
        return String::from("KBBZ");
    }
    // 主诊断在该ADRG主诊断表中, 主手术在该ADRG主手术表中
    let verb_opt: String = adrg_name.to_string() + "_opt";
    let verb_dis: String = adrg_name.to_string() + "_dis";
    return if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && adrg_dis_opt_dict[&verb_opt].contains(&record.main_opt)
    {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组ADRG => 主诊断 + 手术(主手术或其他手术) IB1需要主诊断和两个手术
fn is_both_mdis_opt_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 如果病例没有主手术不符合入组条件
    if record.no_surgery() {
        return String::from("KBBZ");
    }
    // 如果病例没有其他手术也不符合入组条件
    if record.no_other_surgery() {
        return String::from("KBBZ");
    }

    let verb_opt1: String = adrg_name.to_string() + "_op1";
    let verb_opt2: String = adrg_name.to_string() + "_opt2";
    let verb_dis: String = adrg_name.to_string() + "_dis";
    return if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && !(adrg_dis_opt_dict[&verb_opt1].is_disjoint(&record.all_opt))
        && !(adrg_dis_opt_dict[&verb_opt2].is_disjoint(&record.all_opt))
    {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组ADRG => 两个手术
fn is_both_opt_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 如果病例没有其他手术, 不符合入组条件
    if record.no_other_surgery() {
        return String::from("KBBZ");
    }
    let verb_opt1: String = adrg_name.to_string() + "_opt1";
    let verb_opt2: String = adrg_name.to_string() + "_opt2";
    return if !(adrg_dis_opt_dict[&verb_opt1].is_disjoint(&record.all_opt))
        && !(adrg_dis_opt_dict[&verb_opt2].is_disjoint(&record.all_opt))
    {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组的ADRG => 满足主诊断或其他诊断(用于处理PS1\PS2\PS3\PS4)
fn is_dis_adrg(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    let verb_dis: String = adrg_name.to_string() + "_dis";
    if !(adrg_dis_opt_dict[&verb_dis].is_disjoint(&record.all_dis)) {
        return if adrg_name == "PS1" && record.weight < 1500 {
            "PS1".to_string()
        } else if adrg_name == "PS2" && record.weight >= 1500 && record.weight < 1999 {
            "PS2".to_string()
        } else if adrg_name == "PS3" && record.weight >= 1999 && record.weight < 2499 {
            "PS3".to_string()
        } else {
            "PS4".to_string()
        };
    }
    return "KBBZ".to_string();
}

// 特殊ADRG入组, 主诊断+手术表1+手术表2, 或主诊断+手术表1+手术表3+手术表4
fn is_mdis_and_multi_surgery_adrg_one(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 无其他手术则不符合入组条件
    if record.no_other_surgery() {
        return String::from("KBBZ");
    }

    let verb_dis = adrg_name.to_string() + "_dis";
    let verb_opt1 = adrg_name.to_string() + "_opt1";
    let verb_opt2 = adrg_name.to_string() + "_opt2";
    let verb_opt3 = adrg_name.to_string() + "_opt3";
    let verb_opt4 = adrg_name.to_string() + "_opt4";

    return if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && !(adrg_dis_opt_dict[&verb_opt1].is_disjoint(&record.all_opt))
        && !(adrg_dis_opt_dict[&verb_opt2].is_disjoint(&record.all_opt))
    {
        adrg_name
    } else if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && !(adrg_dis_opt_dict[&verb_opt1].is_disjoint(&record.all_opt))
        && !(adrg_dis_opt_dict[&verb_opt3].is_disjoint(&record.all_opt))
        && !(adrg_dis_opt_dict[&verb_opt4].is_disjoint(&record.all_opt))
    {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊ADRG入组, 主诊断+手术表1, 或主诊断+手术表2+手术表3
fn is_mdis_and_multi_surgery_adrg_two(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    let verb_dis: String = adrg_name.to_string() + "_dis";
    let verb_opt1: String = adrg_name.to_string() + "_opt1";
    let verb_opt2: String = adrg_name.to_string() + "_opt2";
    let verb_opt3: String = adrg_name.to_string() + "_opt3";

    return if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && !(adrg_dis_opt_dict[&verb_opt1].is_disjoint(&record.all_opt))
    {
        adrg_name
    } else if adrg_dis_opt_dict[&verb_dis].contains(&record.main_dis)
        && adrg_dis_opt_dict[&verb_opt2].is_disjoint(&record.all_opt)
        && adrg_dis_opt_dict[&verb_opt3].is_disjoint(&record.all_opt)
    {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组的ADRG, 包含全部手术
fn is_all_surgery(record: &DrgCase, all_opt_list: &HashSet<String>, adrg_name: String) -> String {
    // 如果病例无手术则不符合入组条件
    if record.no_surgery() {
        return String::from("KBBZ");
    }

    // 病例手术在所有手术中满足入组条件
    return if !(all_opt_list.is_disjoint(&record.all_opt)) {
        adrg_name
    } else {
        String::from("KBBZ")
    };
}

// 特殊入组的ADRG, 不包含手术
fn is_without_surgery(
    record: &DrgCase,
    all_opt_list: &HashSet<String>,
    adrg_name: String,
) -> String {
    // 病例无手术则直接满足入组条件
    if record.no_surgery() {
        return adrg_name;
    }

    // 病例的手术都不在所有手术列表中则符合入组条件
    return if !(all_opt_list.is_disjoint(&record.all_opt)) {
        String::from("KBBZ")
    } else {
        adrg_name
    };
}

// 特殊入组的ADRG, 没有WB1手术表中的手术
fn is_without_wb1_surgery(
    record: &DrgCase,
    adrg_dis_opt_dict: &HashMap<String, HashSet<String>>,
    adrg_name: String,
) -> String {
    // 病例无手术则不满足入组条件
    if record.no_surgery() {
        return String::from("KBBZ");
    }
    if !(adrg_dis_opt_dict["WB1_opt"].is_disjoint(&record.all_opt)) {
        return adrg_name;
    }
    return String::from("KBBZ");
}

// 综合处理各个ADRG
fn process_adrg(
    record: &DrgCase,
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    all_opt_list: &HashSet<String>,                  // 全部手术列表
    adrg_type_dict: &HashMap<String, Vec<String>>,   // ADRG类型及对应入组类型
    adrg_name: String,
) -> String {
    let pred_adrg = match adrg_type_dict[&adrg_name][2].as_str() {
        // 主手术入组
        "common_opt" => is_common_surgery_adrg(record, adrg_dis_opt, adrg_name),
        // 主诊断入组
        "common_dis" => is_common_diagnosis_adrg(record, adrg_dis_opt, adrg_name),
        // 同时两个手术入组
        "both_opt" => is_both_opt_adrg(record, adrg_dis_opt, adrg_name),
        // 主诊断与主手术入组
        "dis_and_opt" => is_both_mdis_and_mopt_adrg(record, adrg_dis_opt, adrg_name),
        // 主诊断与手术入组
        "main_dis_and_any_opt" => is_both_mdis_opt_adrg(record, adrg_dis_opt, adrg_name),
        // 主诊断与多个手术入组(主诊断+手术表1+手术表2, 主诊断+手术表1+手术表3+手术表4
        "main_dis_and_multi_opt" => {
            is_mdis_and_multi_surgery_adrg_one(record, adrg_dis_opt, adrg_name)
        }
        // 主诊断与多个手术入组(主诊断+手术表1, 主诊断+手术表2+手术表3
        "main_dis_and_multi_opt2" => {
            is_mdis_and_multi_surgery_adrg_two(record, adrg_dis_opt, adrg_name)
        }
        // 两个以上诊断入组(PS1|PS2|PS3|PS4)
        "any_dis" => is_dis_adrg(record, adrg_dis_opt, adrg_name),
        // 所有手术入组
        "all_opt" => is_all_surgery(record, all_opt_list, adrg_name),
        // 无手术入组
        "no_opt" => is_without_surgery(record, all_opt_list, adrg_name),
        // 没有WB1手术入组
        "exclude_wb1_opt" => is_without_wb1_surgery(record, adrg_dis_opt, adrg_name),
        _ => "KBBZ".to_string(),
    };
    return pred_adrg;
}

// 整体分组===========================================================================================
fn which_adrg(
    record: &DrgCase,
    mdc_dis: &HashMap<String, HashSet<String>>, // MDC主诊断表
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    mdcz_main_dis_dict: &HashMap<String, HashSet<String>>, // MDCZ诊断表
    adrg_opt_list: &HashSet<String>,            // 所有手术列表
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG入组类型
    mdc_to_adrg: &HashMap<String, HashSet<String>>, // MDC下的ADRG
) -> (String, String) {
    let mut result_mdc = String::from("KBBZ");
    let mut result_adrg = String::from("KBBZ");
    // 顺序为先期分组 -> 新生儿组 -> 艾滋病组 -> 多发创伤组
    for mdc in [
        "MDCA", "MDCP", "MDCY", "MDCZ", "MDCB", "MDCC", "MDCD", "MDCE", "MDCF", "MDCG", "MDCH",
        "MDCI", "MDCJ", "MDCK", "MDCL", "MDCM", "MDCN", "MDCO", "MDCQ", "MDCR", "MDCS", "MDCT",
        "MDCU", "MDCV", "MDCW", "MDCX",
    ] {
        // 先期分组
        if mdc == "MDCA" {
            result_adrg = is_mdca(
                record,
                adrg_dis_opt,
                adrg_opt_list,
                adrg_type_dict,
                mdc.to_string(),
            );
            if result_adrg != "KBBZ" {
                result_mdc = "MDCA".to_string();
                return (result_adrg, result_mdc);
            }
        }
        // 高优先级MDCP
        else if mdc == "MDCP" {
            result_mdc = is_age_mdc(record, mdc_dis, mdc.to_string());
            if result_mdc != "KBBZ" {
                for adrg in &mdc_to_adrg[mdc] {
                    result_adrg = process_adrg(
                        record,
                        adrg_dis_opt,
                        adrg_opt_list,
                        adrg_type_dict,
                        adrg.to_string()
                    );
                    if result_adrg != "KBBZ" {
                        return (result_adrg, result_mdc);
                    }
                }
            }
        }
        // 高优先级MDCZ
        else if mdc == "MDCZ" {
            result_mdc = is_mdcz(record, mdcz_main_dis_dict, mdc.to_string());
            if result_mdc != "KBBZ" {
                for adrg in &mdc_to_adrg[mdc] {
                    result_adrg = process_adrg(
                        record,
                        adrg_dis_opt,
                        adrg_opt_list,
                        adrg_type_dict,
                        adrg.to_string()
                    );
                    if result_adrg != "KBBZ" {
                        return (result_adrg, result_mdc);
                    }
                }
            }
        }
        // 需要考虑性别的MDC大类
        else if mdc == "MDCM" || mdc == "MDCN" {
            result_mdc = is_sex_mdc(record, mdc_dis, mdc.to_string());
            if result_mdc != "KBBZ" {
                for adrg in &mdc_to_adrg[mdc] {
                    result_adrg = process_adrg(
                        record,
                        adrg_dis_opt,
                        adrg_opt_list,
                        adrg_type_dict,
                        adrg.to_string()
                    );
                    if result_adrg != "KBBZ" {
                        return (result_adrg, result_mdc);
                    }
                }
            }
        }
        // 普通凭借主诊断入组的MDC大类
        else {
            result_mdc = is_common_mdc(record, mdc_dis, mdc.to_string());
            if result_mdc != "KBBZ" {
                for adrg in &mdc_to_adrg[mdc] {
                    result_adrg = process_adrg(
                        record,
                        adrg_dis_opt,
                        adrg_opt_list,
                        adrg_type_dict,
                        adrg.to_string()
                    );
                    if result_adrg != "KBBZ" {
                        return (result_adrg, result_mdc);
                    }
                }
            }
        }
    }

    return (result_adrg, result_mdc);
}

// 判断是否为QY病例
fn is_qy(
    record: &DrgCase,                              // 病例结构
    adrg_pred: String,                             // 已经入的ADRG组
    mdc_pred: String,                              // 已经进入的MDC大类
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG入组类型字典
    all_opt_list: &HashSet<String>,                // 所有手术列表
) -> (String, String) {
    if mdc_pred != "KBBZ" {
        // 包含全部手术的, 不会出现QY
        if adrg_pred == "YC1".to_string()
            || adrg_pred == "SB1".to_string()
            || adrg_pred == "XJ1".to_string()
            || adrg_pred == "TB1".to_string()
        {
            return (adrg_pred, mdc_pred);
        }
        // 无手术的病例也不会出现QY
        if record.no_surgery() {
            return (adrg_pred, mdc_pred);
        }
        // 初分组为内科ADRG且有有效手术的病例为QY病例
        if (adrg_type_dict[&adrg_pred][0] == "内科") && !(all_opt_list.is_disjoint(&record.all_opt))
        {
            let new_adrg_pred = match mdc_pred.as_str() {
                "MDCA" => "AQY",
                "MDCB" => "BQY",
                "MDCC" => "CQY",
                "MDCD" => "SQY",
                "MDCE" => "EQY",
                "MDCF" => "FQY",
                "MDCG" => "GQY",
                "MDCH" => "HQY",
                "MDCI" => "IQY",
                "MDCJ" => "JQY",
                "MDCK" => "KQY",
                "MDCL" => "LQY",
                "MDCM" => "MQY",
                "MDCN" => "NQY",
                "MDCO" => "OQY",
                "MDCP" => "PQY",
                "MDCQ" => "QQY",
                "MDCR" => "RQY",
                "MDCU" => "UQY",
                "MDCV" => "VQY",
                "MDCW" => "WQY",
                "MDCZ" => "ZQY",
                _ => "KBBZ",
            };
            return (new_adrg_pred.to_string(), mdc_pred);
        }
    }
    return (adrg_pred.to_string(), mdc_pred);
}

// 判断CCMCC
fn cc_mcc(
    record: &DrgCase,                              // 病例结构
    adrg_pred: String,                             // 已经入的ADRG组
    mdc_pred: String,                              // 已经进入的MDC大类
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG类型
    exclude_dict: &HashMap<String, String>,        // 排除表
    cc_mcc_dict: &HashMap<String, Vec<String>>,    // CCMCC表
) -> String {
    let mut complication_list: Vec<String> = Vec::new();
    let mut exclude_pos = "";
    let mut complication = "";

    // 如果该ADRG没有并发症细分, 则并发症类型为9
    if adrg_type_dict[&adrg_pred][1].as_str() == "未细分" {
        return "9".to_string();
    }

    // 如果无其他诊断, 则病例无并发症
    if record.no_other_diagnosis() {
        return "5".to_string();
    }

    let default_dict_val_vec = vec!["".to_string(), "".to_string()];
    let default_dict_val_str = "".to_string();
    // let default_dict_val_str = "";
    // 有其他诊断的情况下, 逐一检查是否为CC或MCC, 是否被排除
    for _d in &record.other_dis {
        let temp = cc_mcc_dict.get(_d).unwrap_or(&default_dict_val_vec);
        exclude_pos = &temp[0];
        complication = &temp[1];
        // 有严重或一般并发症
        if complication != "" {
            if  exclude_pos != "无" {
                // 是否被排除
                let is_exclude = exclude_dict
                    .get(&record.main_dis)
                    .unwrap_or(&default_dict_val_str);
                if is_exclude != exclude_pos {
                    complication_list.push(complication.to_string());
                }
            } else {
                // 没有排除表
                complication_list.push(complication.to_string());
            }

        }
    }
    if adrg_type_dict[&adrg_pred][1].as_str() == "1合并3" {
        if complication_list.len() == 0 {
            return "5".to_string();
        } else {
            return "3".to_string();
        }
    } else if adrg_type_dict[&adrg_pred][1].as_str() == "3合并5" {
        if complication_list.contains(&"MCC".to_string()) {
            return "1".to_string();
        } else {
            return "5".to_string();
        }
    } else {
        if complication.len() == 0 {
            return "5".to_string();
        } else if complication_list.contains(&"MCC".to_string()) {
            return "1".to_string();
        } else {
            return "3".to_string();
        }
    }
}

fn which_drg(
    record: &DrgCase,
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG类型
    exclude_dict: &HashMap<String, String>,        // 排除表
    cc_mcc_dict: &HashMap<String, Vec<String>>,    // CCMCC表
    mdc_dis: &HashMap<String, HashSet<String>>,    // MDC主诊断表
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    mdcz_main_dis_dict: &HashMap<String, HashSet<String>>, // MDCZ诊断表
    adrg_opt_list: &HashSet<String>,               // 所有手术列表
    mdc_to_adrg: &HashMap<String, HashSet<String>>, // MDC下的ADRG
) -> String {
    // 判断进入的MDC
    let (mut adrg, mut mdc) = which_adrg(
        record,
        mdc_dis,
        adrg_dis_opt,
        mdcz_main_dis_dict,
        adrg_opt_list,
        adrg_type_dict,
        mdc_to_adrg,
    );

    // 判断是否为QY
    (adrg, mdc) = is_qy(
        record,
        adrg.clone(),
        mdc.clone(),
        adrg_type_dict,
        adrg_opt_list,
    );
    // 判断CCMCC标志
    if adrg != "KBBZ" && &adrg[1..3] != "QY" {
        let ccmcc_lab = cc_mcc(
            record,
            adrg.clone(),
            mdc.clone(),
            adrg_type_dict,
            exclude_dict,
            cc_mcc_dict,
        );
        return adrg.to_string() + ccmcc_lab.as_str();
    } else {
        return adrg.to_string();
    }
}

// 用户输入======================================================================
// 读取用户输入的文本向量
fn read_vec_from_terminal() -> Vec<String> {
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    if input.trim().is_empty() {
        Vec::new()
    } else {
        input
            .trim()
            .split(',')
            .map(|s| s.trim().to_string())
            .collect()
    }
}
// 读取用户输入的整数
fn read_int_from_terminal() -> i32 {
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    input.trim().parse::<i32>().unwrap()
}

// 读取用户输入的浮点数
fn read_float_from_terminal() -> f64 {
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    input.trim().parse::<f64>().unwrap()
}

// 读取单个字符串
fn read_str_from_terminal() -> String {
    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();
    input.trim().to_string()
}

// 接收用户在命令行的输入, 初始化结构体
fn create_drg_case_from_terminal() -> DrgCase {
    println!("Enter admission number: ");
    let admission_number = read_str_from_terminal();

    println!("Enter principal diagnosis: ");
    let principal_diagnosis = read_str_from_terminal();

    println!("Enter principal operation: ");
    let principal_operation = read_str_from_terminal();

    println!("Enter other diagnosis: ");
    let other_diagnosis = read_vec_from_terminal();

    println!("Enter other operation: ");
    let other_operation = read_vec_from_terminal();

    println!("Enter gender: ");
    let gender = read_int_from_terminal();

    println!("Enter age: ");
    let age = read_float_from_terminal();

    println!("Enter weight: ");
    let mass = read_int_from_terminal();

    DrgCase::new(
        admission_number,
        principal_diagnosis,
        principal_operation,
        other_diagnosis,
        other_operation,
        gender,
        age,
        mass,
    )
}

// 读取表格文件✔
fn from_csv_file(file_path: &str) -> Result<DataFrame, Box<dyn Error>> {
    let res = CsvReader::from_path(file_path)?
        .infer_schema(None)
        .has_header(true)
        .finish()?;
    Ok(res)
}

// 转换ICD编码, 在ICD10中涉及到x的只有小写
fn icd_transform(icd: String) -> String {
    icd.chars().map(|c| {
        if c == 'x' {
            c
        } else {
            c.to_uppercase().collect::<String>().chars().next().unwrap()
        }
    }).collect::<String>()
}

// 合并表中的各个其他手术列的手术编码为一个向量✔
fn concat_icd9_code(df: &DataFrame, idx: usize) -> Result<Vec<String>, Box<dyn Error>> {
    let mut row_val: Vec<String> = Vec::new();
    for colname in ["其他手术编码1","其他手术编码2","其他手术编码3","其他手术编码4","其他手术编码5","其他手术编码6","其他手术编码7",
                    "其他手术编码8","其他手术编码9","其他手术编码10","其他手术编码11","其他手术编码12","其他手术编码13","其他手术编码14",
                    "其他手术编码15","其他手术编码16"]
    {
        if df.column(colname)?.get(idx).unwrap() != AnyValue::Null {
            let res = df.column(colname)?.get(idx)?;
            row_val.push(res.to_string());
        }
    }
    let modify_vec: Vec<String> = row_val.into_iter().map(|s| s.trim_matches('\"').to_string()).collect();
    return Ok(modify_vec)
}

// 合并表格中的各个其他诊断编码为一个向量(需要转换大小写)✔
fn concat_icd10_code(df: &DataFrame, idx: usize) -> Result<Vec<String>, Box<dyn Error>> {
    let mut row_val: Vec<String> = Vec::new();
    for colname in ["其他诊断编码1","其他诊断编码2","其他诊断编码3","其他诊断编码4","其他诊断编码5","其他诊断编码6","其他诊断编码7",
        "其他诊断编码8","其他诊断编码9","其他诊断编码10","其他诊断编码11","其他诊断编码12","其他诊断编码13","其他诊断编码14",
        "其他诊断编码15","其他诊断编码16"]
    {
        if df.column(colname)?.get(idx).unwrap() != AnyValue::Null {
            let res = df.column(colname)?.get(idx)?;
            row_val.push(res.to_string());
        }
    }
    let modify_vec: Vec<String> = row_val.into_iter().map(|s| icd_transform(s.trim_matches('\"').to_string())).collect();
    return Ok(modify_vec)
}

// 从表格数据构造出DRG病例结构
fn construct_drg_case(df: &DataFrame) -> Result<Vec<DrgCase>, Box<dyn Error>> {
    let df_size = df.shape();
    let mut my_vec: Vec<DrgCase> = Vec::new();
    for i in 0..df_size.0 as usize {
        let temp_main_opt = df.column("主手术编码")?.get(i).unwrap();
        let mut drg_case = DrgCase::new(
            df.column("结算流水号")?.get(i)?.to_string(),
            icd_transform(df.column("主诊断编码")?.get(i)?.to_string()),
            match temp_main_opt {
                AnyValue::Null => "".to_string(),
                _ => temp_main_opt.to_string()
            },
            concat_icd10_code(df, i)?,
            concat_icd9_code(df, i)?,
            df.column("性别")?.i32()?.get(i).unwrap(),
            df.column("年龄")?.f64()?.get(i).unwrap(),
            df.column("体重")?.i32()?.get(i).unwrap()
        );
        let _ = drg_case.concat_dis();   // 将其他诊断与主诊断合并在一起
        let _ = &drg_case.concat_opt();  // 将其他手术与主手术合并在一起
        my_vec.push(drg_case);
    }
    return Ok(my_vec)
}

// 批量对表格数据进行DRG分组
fn batch_drg_group(
    df: &DataFrame,
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG类型
    exclude_dict: &HashMap<String, String>,        // 排除表
    cc_mcc_dict: &HashMap<String, Vec<String>>,    // CCMCC表
    mdc_dis: &HashMap<String, HashSet<String>>,    // MDC主诊断表
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    mdcz_main_dis_dict: &HashMap<String, HashSet<String>>, // MDCZ诊断表
    adrg_opt_list: &HashSet<String>,               // 所有手术列表
    mdc_to_adrg: &HashMap<String, HashSet<String>>, // MDC下的ADRG
) -> Vec<String> {
    let mut pred_drg_list: Vec<String> = Vec::new();
    let drg_case = construct_drg_case(&df).unwrap();
    for case in drg_case {
        let drg_pred = which_drg(
            &case,
            adrg_type_dict,
            exclude_dict,
            cc_mcc_dict,
            mdc_dis,
            adrg_dis_opt,
            mdcz_main_dis_dict,
            adrg_opt_list,
            mdc_to_adrg
        );
        pred_drg_list.push(drg_pred);
    }
    return pred_drg_list
}

// 对表格数据进行DRG分组并导出原表格及分组结果
fn drg_group_and_export(
    in_path: &str,
    out_path: &str,
    adrg_type_dict: &HashMap<String, Vec<String>>, // ADRG类型
    exclude_dict: &HashMap<String, String>,        // 排除表
    cc_mcc_dict: &HashMap<String, Vec<String>>,    // CCMCC表
    mdc_dis: &HashMap<String, HashSet<String>>,    // MDC主诊断表
    adrg_dis_opt: &HashMap<String, HashSet<String>>, // ADRG诊断手术表
    mdcz_main_dis_dict: &HashMap<String, HashSet<String>>, // MDCZ诊断表
    adrg_opt_list: &HashSet<String>,               // 所有手术列表
    mdc_to_adrg: &HashMap<String, HashSet<String>>, // MDC下的ADRG
) -> Result<(), Box<dyn Error>> {
    // 读取CSV表格文件
    let mut df = from_csv_file(in_path).unwrap();
    // 进行DRG分组
    let drg_pred_list = batch_drg_group(&df,
                                        adrg_type_dict,
                                        exclude_dict,
                                        cc_mcc_dict,
                                        mdc_dis,
                                        adrg_dis_opt,
                                        mdcz_main_dis_dict,
                                        adrg_opt_list,
                                        mdc_to_adrg);
    // 创建一个Series序列准备添加到数据表中
    let new_col = Series::new("clear_code", drg_pred_list);
    // 定义需要添加的列的数据类型
    let mut my_schema = Schema::new();
    my_schema.with_column(String::from("clear_code"), DataType::UInt8);
    // 向表中添加列
    &df._add_columns(vec![new_col], &my_schema)?;

    // 将表格数据以CSV格式写入本地
    let export_file = File::create(&out_path)?;
    CsvWriter::new(export_file)
        .has_header(true)
        .with_delimiter(b',')
        .finish(&mut df)?;

    Ok(())
}

fn main() {
    // 读取分组方案===================================================================================
    // 读取MDC主诊断表
    let mdc_main_dis_dict = read_file_as_str_to_set("data\\MDC_main_dis.json").unwrap();
    // 读取MDCZ主诊断表
    let mdcz_main_dis_dict = read_file_as_str_to_set("data\\MDCZ_main_dis_list.json").unwrap();
    // 读取ADRG诊断手术表
    let adrg_dis_opt_dict = read_file_as_str_to_set("data\\adrg_dis_opt.json").unwrap();
    // 读取ADRG类型表
    let adrg_drg_dict = read_file_as_str_to_tuple("data\\adrg_type_dict.json").unwrap();
    // 读取MDC对应ADRG表
    let mdc_to_adrg_dict = read_file_as_str_to_set("data\\mdc_map_adrg.json").unwrap();
    // 读取CCMCC表
    let cc_mcc_dict = read_file_as_str_to_tuple("data\\cc_mcc_dict.json").unwrap();
    // 读取排除表
    let exclusive_dict = read_file_as_str_to_str("data\\exclusive_dict.json").unwrap();
    // 读取所有手术列表
    let all_adrg_opt_list = read_icd9_to_vec("data\\all_opt_sheet.txt").unwrap();

    // 构造病例结构===================================================================================
    // let mut test_case = DrgCase::new(
    //     String::from("450100G0000013706219"),
    //     String::from("J20.900"),
    //     String::from("93.3500x004"),
    //     vec!["E87.102".to_string(), "E87.803".to_string()],
    //     vec![],
    //     1,
    //     29.0,
    //     2789,
    // );
    // let _ = &test_case.concat_opt();
    // let _ = &test_case.concat_dis();

    // let drg_pred = which_drg(
    //     &test_case,
    //     &adrg_drg_dict,
    //     &exclusive_dict,
    //     &cc_mcc_dict,
    //     &mdc_main_dis_dict,
    //     &adrg_dis_opt_dict,
    //     &mdcz_main_dis_dict,
    //     &all_adrg_opt_list,
    //     &mdc_to_adrg_dict,
    // );

    // 选择模式(单个病例分组或输入表格批量分组)
    let mut selected_mood = "single".to_string();
    println!("Please choose mood: [single] or [batch]");
    selected_mood = read_str_from_terminal();

    if selected_mood == "single" {
        // 单个分组
        let mut go_on = "yes".to_string();

        while go_on == "yes" {
            // println!("please enter infomation about case: ");
            let mut this_case = create_drg_case_from_terminal();
            let _ = &this_case.concat_dis();
            let _ = &this_case.concat_opt();
            let this_drg_pred = which_drg(
                &this_case,
                &adrg_drg_dict,
                &exclusive_dict,
                &cc_mcc_dict,
                &mdc_main_dis_dict,
                &adrg_dis_opt_dict,
                &mdcz_main_dis_dict,
                &all_adrg_opt_list,
                &mdc_to_adrg_dict,
            );
            println!("{:?}", &this_case.all_opt);
            println!("{:?}", &this_case.all_dis);
            println!("Drg code of this case is {}", this_drg_pred);
            println!("Do you want continue drg group? enter yes to continue, otherwise enter quit");
            go_on = read_str_from_terminal();
        }
    } else {
        // 批量表格分组
        println!("please enter import file path: ");
        let in_file_path= read_str_from_terminal();
        println!("please enter export file path: ");
        let out_file_path = read_str_from_terminal();
        drg_group_and_export(
            &in_file_path,
            &out_file_path,
            &adrg_drg_dict,
            &exclusive_dict,
            &cc_mcc_dict,
            &mdc_main_dis_dict,
            &adrg_dis_opt_dict,
            &mdcz_main_dis_dict,
            &all_adrg_opt_list,
            &mdc_to_adrg_dict
        ).expect("drg group fail please check if there are any wrong in dataset");
    }

}
