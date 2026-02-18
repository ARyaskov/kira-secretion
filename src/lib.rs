pub mod aggregate;
pub mod cli;
pub mod expr;
pub mod input;
pub mod model;
pub mod panels;
pub mod pipeline;
pub mod report;
pub mod simd;

pub mod prelude {
    pub use crate::input::detect::TenXFormat;
    pub use crate::pipeline::stage1_load::DatasetCtx;
}
