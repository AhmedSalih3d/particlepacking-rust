use crate::itertools::{Itertools, Product};

use smallvec::SmallVec;
use std::iter::Iterator;

use nalgebra::{distance_squared, Point3};

const DIM: usize = 3;

type CellIdx = usize;

type EntityIdx = usize;

pub struct GridPartition {
    // squared search radius
    squared_radius: f32,
    // number of subdivisions for each dimension
    pub nb_subdivisons: [usize; DIM],
    // length of the cells for each dimension (slightly larger than the search radius)
    pub cell_width: [f32; DIM],
    // total number of cells
    pub nb_cells: usize,
    // adjascent_cells[cell_idx] = indexes of all cells adjascent to 'cell_idx'
    // adjascent_cells.len() = #cells
    pub adjascent_cells: Vec<SmallVec<[CellIdx; 14]>>,
    // head[cell_idx] = position of the first element of the cell in self.successor
    // head.len() = #cells.
    head: Vec<Option<EntityIdx>>,
    // successor[entity_idx] = idx of the sucessor of entity 'entity_idx' in the cell
    // successor.len() = #entities.
    successor: Vec<Option<EntityIdx>>,
}

impl GridPartition {
    // Create a new grid, provided a search radius
    // Note: the points are assumed o belong to 
    // [0, dimensions[0]] x [0, dimensions[1]] x [0, dimensions[2]] 
    pub fn new(dimensions: [f32; DIM], radius: f32) -> Self {
        let squared_radius = radius.powi(2);
        let nb_subdivisons = [
            (dimensions[0] / radius) as usize,
            (dimensions[1] / radius) as usize,
            (dimensions[2] / radius) as usize,
        ];
        let cell_width = [
            dimensions[0] / nb_subdivisons[0] as f32,
            dimensions[1] / nb_subdivisons[1] as f32,
            dimensions[2] / nb_subdivisons[2] as f32,
        ];
        let nb_cells = (nb_subdivisons[0] * nb_subdivisons[1] * nb_subdivisons[2]) as usize;
        let adjascent_cells = vec![SmallVec::<[CellIdx; 14]>::new(); nb_cells];
        let head = vec![None; nb_cells];
        let successor = Vec::new();
        let mut grid = Self {
            squared_radius,
            nb_subdivisons,
            cell_width,
            nb_cells,
            adjascent_cells,
            head,
            successor,
        };
        grid.build_adjascent_cells();
        grid
    }

    // Precompute the 'adjascent_cells' vector
    fn build_adjascent_cells(&mut self) {
        // for every cell
        for idx_flat in 0..self.nb_cells {
            let idx_spatial = self.unflatten_index(idx_flat);
            // for all adjascent cells
            for &dx in &[-1, 0, 1] {
                for &dy in &[-1, 0, 1] {
                    for &dz in &[-1, 0, 1] {
                        // if the adjascent cell is not out of bounds
                        if (idx_spatial.x as i64 + dx >= 0)
                            && (idx_spatial.x as i64 + dx < self.nb_subdivisons[0] as i64)
                            && (idx_spatial.y as i64 + dy >= 0)
                            && (idx_spatial.y as i64 + dy < self.nb_subdivisons[1] as i64)
                            && (idx_spatial.z as i64 + dz >= 0)
                            && (idx_spatial.z as i64 + dz < self.nb_subdivisons[2] as i64)
                        {
                            let other_idx = Point3::new(
                                (idx_spatial.x as i64 + dx) as usize,
                                (idx_spatial.y as i64 + dy) as usize,
                                (idx_spatial.z as i64 + dz) as usize,
                            );
                            let other_idx_flat = self.flatten_index(other_idx);
                            // makes sure to avoid duplicate neighbors: two cells (i, j)
                            // are adjascent only if i <=j
                            if other_idx_flat <= idx_flat {
                                // add the adjascent cell to the vector
                                self.adjascent_cells[idx_flat].push(other_idx_flat);
                            }
                        }
                    }
                }
            }
        }
    }

    // given the sptial coordinates (i, j, k)  of a cell, compute its integer index.
    fn flatten_index(&self, idx_spatial: Point3<usize>) -> usize {
        idx_spatial.x * self.nb_subdivisons[1] * self.nb_subdivisons[2]
            + idx_spatial.y * self.nb_subdivisons[2]
            + idx_spatial.z
    }

    // given the integer index of a cell, compute its (i, j, k) coordinates
    fn unflatten_index(&self, idx_flat: CellIdx) -> Point3<usize> {
        Point3::new(
            idx_flat / (self.nb_subdivisons[1] * self.nb_subdivisons[2]),
            (idx_flat / self.nb_subdivisons[2]) % self.nb_subdivisons[1],
            idx_flat % self.nb_subdivisons[2],
        )
    }

    // Return an iterator over the entities contained in a cell
    pub fn iter_cell(&self, idx_flat: CellIdx) -> CellIterator {
        CellIterator {
            current_entity: self.head[idx_flat],
            successor: &self.successor,
        }
    }

    // Return an iterator over the entities contained in cells adjascent to cell 'idx_flat'
    pub fn iter_adjascent_cells(&self, idx_flat: CellIdx) -> MultiCellIterator {
        let adjascent_cells = &self.adjascent_cells[idx_flat];
        let current_cell = Some(self.iter_cell(adjascent_cells[0]));
        let remaining_cells: SmallVec<[CellIterator; 14]> = adjascent_cells[1..]
            .iter()
            .map(|&i| self.iter_cell(i))
            .collect();
        MultiCellIterator {
            current_cell,
            remaining_cells,
        }
    }

    // Return an iterator over all couples the couples (i, j), where j is in cell 'idx_flat',
    // and i is in an cell adjascent to 'idx_cell'
    fn iter_cell_couples(&self, idx_flat: CellIdx) -> Product<MultiCellIterator, CellIterator> {
        self.iter_adjascent_cells(idx_flat)
            .cartesian_product(self.iter_cell(idx_flat))
    }

    // Return an iterator over all couples (i, j), where i and j are within the (squared) distance self.squared_radius
    pub fn iter_neighbors<'a>(&'a self, points: &'a [Point3<f32>]) -> NeighborIterator {
        NeighborIterator {
            grid: self,
            current_cell: 0,
            cartesian_iterator: self.iter_cell_couples(0),
            points,
        }
    }

    // Reset the head and the successor vector
    fn empty_grid(&mut self, nb_points: usize) {
        // Reset the 'head' vector
        for x in &mut self.head {
            *x = None;
        }
        // Resize the 'successor' vector if needed
        let additional_capacity = nb_points - self.successor.capacity();
        if additional_capacity > 0 {
            self.successor.reserve_exact(additional_capacity)
        }
        // Fill the successor vector with None values.
        // the values already present are not replaced,
        // since they will be overwritten in self.fill
        let additional_points = nb_points - self.successor.len();
        if additional_points > 0 {
            for _ in 0..additional_points {
                self.successor.push(None);
            }
        }
    }

    // fill the grid with a new set of points.
    fn fill_grid(&mut self, points: &[Point3<f32>]) {
        self.empty_grid(points.len());
        // for each point
        for (i, point) in points.iter().enumerate() {
            let idx_spatial = Point3::new(
                (point.x / self.cell_width[0]) as usize,
                (point.y / self.cell_width[1]) as usize,
                (point.z / self.cell_width[2]) as usize,
            );
            let idx_flat = self.flatten_index(idx_spatial);
            // define the successor of point as the previous first element of the cell
            self.successor[i] = self.head[idx_flat];
            // define the first element of the cell as the current point
            self.head[idx_flat] = Some(i);
        }
    }

    // Update the grid with a new set of points,
    // and return an iterator over all the neighbor indexes
    pub fn query_neighbors<'a>(&'a mut self, points: &'a [Point3<f32>]) -> NeighborIterator {
        self.fill_grid(points);
        self.iter_neighbors(points)
    }
}

#[derive(Clone)]
// Iterator over all the entities in a cell
pub struct CellIterator<'a> {
    current_entity: Option<EntityIdx>,
    successor: &'a Vec<Option<EntityIdx>>,
}

impl<'a> Iterator for CellIterator<'a> {
    type Item = EntityIdx;

    fn next(&mut self) -> Option<EntityIdx> {
        // if their is a current entity
        if let Some(i) = self.current_entity {
            // get the next entity
            self.current_entity = self.successor[i];
            // return the current entity
            Some(i)
        // else: we finished iterating the cell
        } else {
            None
        }
    }
}

// Iterator over the entities of several cells
pub struct MultiCellIterator<'a> {
    // Would it be better to use an iterator for 'remaining_cells' instead of an array ? How to do it ?
    current_cell: Option<CellIterator<'a>>,
    remaining_cells: SmallVec<[CellIterator<'a>; 14]>,
}

impl<'a> Iterator for MultiCellIterator<'a> {
    type Item = EntityIdx;

    fn next(&mut self) -> Option<EntityIdx> {
        // if their is a current cell to iterate over
        if let Some(cell) = &mut self.current_cell {
            // if the current cell is not empty
            if let Some(i) = &mut cell.next() {
                // pop the element of the cell
                Some(*i)
            // else: we finished iterating the current cell
            } else {
                // go to the next cell
                self.current_cell = self.remaining_cells.pop();
                // keep iterating
                self.next()
            }
        // else: we have finished visiting all cells
        } else {
            None
        }
    }
}

// Iterator over all neighbors of the Grid
pub struct NeighborIterator<'a> {
    grid: &'a GridPartition,
    current_cell: CellIdx,
    cartesian_iterator: itertools::Product<MultiCellIterator<'a>, CellIterator<'a>>,
    points: &'a [Point3<f32>],
}

impl Iterator for NeighborIterator<'_> {
    type Item = (CellIdx, EntityIdx, f32);

    fn next(&mut self) -> Option<(EntityIdx, EntityIdx, f32)> {
        // if the cartesin iterator is not empty
        if let Some((idx, idx_other)) = self.cartesian_iterator.next() {
            if idx < idx_other {
                // if the points are close enought, yield the neighbors and their distance
                let squared_distance =
                    distance_squared(&self.points[idx], &self.points[idx_other]) as f32;
                if squared_distance < self.grid.squared_radius {
                    Some((idx, idx_other, squared_distance.powf(0.5)))
                } else {
                    self.next()
                }
            // otherwise, keep iterating
            } else {
                self.next()
            }
        // else: the cartesian iterator is empty
        } else {
            // go to the next cell
            self.current_cell += 1;
            // if the current cell is not the last one
            if self.current_cell < self.grid.nb_cells {
                // compute the cartesian iterator of the next cell
                self.cartesian_iterator = self.grid.iter_cell_couples(self.current_cell);
                // keep iterating
                self.next()
            // else:
            } else {
                // the cartesian iterator is empty this is the last cell
                // we have finised iterating
                None
            }
        }
    }
}
