//! Experimental construction-topology kernel, separate from the legacy CLI.
pub mod reconstruction;
// Share input types and parser; v2 does not call the legacy reconstructor.
pub mod input;
pub mod parsers;
