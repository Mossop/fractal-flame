use std::fmt::Display;

#[cfg(feature = "logging")]
pub type RunState = OrderedKV;
#[cfg(not(feature = "logging"))]
pub type RunState = bool;

pub trait LogRepr {
    fn repr(self) -> String;
}

impl LogRepr for String {
    fn repr(self) -> String {
        self
    }
}

impl LogRepr for &str {
    fn repr(self) -> String {
        self.to_string()
    }
}

impl LogRepr for f64 {
    fn repr(self) -> String {
        format!("{:.20}", self)
    }
}

impl LogRepr for usize {
    fn repr(self) -> String {
        format!("{}", self)
    }
}

impl LogRepr for u32 {
    fn repr(self) -> String {
        format!("{}", self)
    }
}

impl LogRepr for i32 {
    fn repr(self) -> String {
        format!("{}", self)
    }
}

#[derive(Default, Clone)]
pub struct OrderedKV {
    pairs: std::collections::HashMap<String, String>,
    order: Vec<String>,
}

impl OrderedKV {
    pub fn insert<V: LogRepr>(&mut self, key: String, value: V) {
        if !self.pairs.contains_key(&key) {
            self.order.push(key.clone());
        }

        self.pairs.insert(key, value.repr());
    }
}

impl Display for OrderedKV {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut elements: Vec<String> = Vec::new();
        for key in &self.order {
            elements.push(format!("{}={}", key, self.pairs.get(key).unwrap()));
        }

        f.pad(&elements.join(", "))
    }
}

#[cfg(feature = "logging")]
macro_rules! state {
    ({
        $(
            $field:ident: $value:expr,
        )*
    }) => {{
        let mut state = crate::logging::RunState::default();
        $(
            state.insert(stringify!($field).to_string(), $value);
        )*
        state
    }};
    ($state:expr, {
        $(
            $field:ident: $value:expr,
        )*
    }) => {{
        let mut state = $state.clone();
        $(
            state.insert(stringify!($field).to_string(), $value);
        )*
        state
    }};
}
#[cfg(not(feature = "logging"))]
macro_rules! state {
    ({
        $(
            $field:ident: $value:expr,
        )*
    }) => {
        false
    };
    ($state:expr, {
        $(
            $field:ident: $value:expr,
        )*
    }) => {
        false
    };
}

#[cfg(feature = "logging")]
macro_rules! runlog {
    ($message:expr, $state:expr) => {
        println!("{}: {}", $message, $state);
    };
    ($state:expr) => {
        println!("{}", $state);
    };
}
#[cfg(not(feature = "logging"))]
macro_rules! runlog {
    ($message:expr, $state:expr) => {};
    ($state:expr) => {};
}

pub(crate) use {runlog, state};
