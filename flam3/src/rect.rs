use std::ops::{Index, IndexMut};

pub struct Rect<T> {
    width: usize,
    buffer: Vec<T>,
}

impl<T> Rect<T> {
    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.buffer.len() / self.width
    }

    pub fn buffer(&self) -> &[T] {
        &self.buffer
    }

    pub fn buffer_mut(&mut self) -> &mut [T] {
        &mut self.buffer
    }

    pub fn buffer_index(&self, x: usize, y: usize) -> usize {
        self.width * y + x
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.buffer.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.buffer.iter_mut()
    }
}

impl<T: Default + Clone> Rect<T> {
    pub fn default(width: usize, height: usize) -> Rect<T> {
        Rect {
            width,
            buffer: vec![Default::default(); width * height],
        }
    }
}

impl<T> Index<(usize, usize)> for Rect<T> {
    type Output = T;

    fn index(&self, (x, y): (usize, usize)) -> &T {
        &(self.buffer()[self.buffer_index(x, y)])
    }
}

impl<T> IndexMut<(usize, usize)> for Rect<T> {
    fn index_mut(&mut self, (x, y): (usize, usize)) -> &mut T {
        let index = self.buffer_index(x, y);
        &mut (self.buffer_mut()[index])
    }
}

impl<T> Index<usize> for Rect<T> {
    type Output = T;

    fn index(&self, idx: usize) -> &T {
        self.buffer.index(idx)
    }
}

impl<T> IndexMut<usize> for Rect<T> {
    fn index_mut(&mut self, idx: usize) -> &mut T {
        self.buffer.index_mut(idx)
    }
}
