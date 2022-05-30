#![allow(dead_code)]

pub type RunState = OrderedKV;

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
    #[cfg(feature = "logging")]
    pairs: std::collections::HashMap<String, String>,
    #[cfg(feature = "logging")]
    order: Vec<String>,
}

impl OrderedKV {
    #[cfg(feature = "logging")]
    pub fn insert<V: LogRepr>(&mut self, key: String, value: V) {
        if !self.pairs.contains_key(&key) {
            self.order.push(key.clone());
        }

        self.pairs.insert(key, value.repr());
    }

    #[cfg(not(feature = "logging"))]
    pub fn insert<V: LogRepr>(&mut self, _key: String, _value: V) {}
}

#[cfg(feature = "logging")]
impl std::fmt::Display for OrderedKV {
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
        crate::logging::RunState::default()
    };
    ($state:expr, {
        $(
            $field:ident: $value:expr,
        )*
    }) => {
        $state.clone()
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
