use crate::basis::Basis;
use crate::graph::Graph;
use crate::types::{BasisMatrix, CoordinateVector, OffsetVector};
use crate::unitcell::Unitcell;
use roxmltree::{Document, Node};
use std::fmt::{self, Write as _};
use std::fs;
use std::path::Path;

#[derive(Debug)]
pub enum XmlError {
  Parse(roxmltree::Error),
  Io(std::io::Error),
  MissingAttribute(&'static str),
  InvalidFormat(&'static str),
  DimensionMismatch(&'static str),
  UnknownEntry(String),
}

impl fmt::Display for XmlError {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      XmlError::Parse(err) => write!(f, "XML parse error: {err}"),
      XmlError::Io(err) => write!(f, "I/O error: {err}"),
      XmlError::MissingAttribute(attr) => write!(f, "missing XML attribute: {attr}"),
      XmlError::InvalidFormat(msg) => write!(f, "invalid XML format: {msg}"),
      XmlError::DimensionMismatch(msg) => write!(f, "dimension mismatch: {msg}"),
      XmlError::UnknownEntry(name) => write!(f, "lattice entry not found: {name}"),
    }
  }
}

impl std::error::Error for XmlError {}

impl From<roxmltree::Error> for XmlError {
  fn from(value: roxmltree::Error) -> Self {
    XmlError::Parse(value)
  }
}

impl From<std::io::Error> for XmlError {
  fn from(value: std::io::Error) -> Self {
    XmlError::Io(value)
  }
}

pub type Result<T> = std::result::Result<T, XmlError>;

pub fn read_basis_from_str(xml: &str, name: &str) -> Result<Basis> {
  let doc = Document::parse(xml)?;
  let entry = find_entry(&doc, "LATTICE", name)?.ok_or_else(|| XmlError::UnknownEntry(name.to_string()))?;
  parse_basis(entry)
}

pub fn read_unitcell_from_str(xml: &str, name: &str) -> Result<Unitcell> {
  let doc = Document::parse(xml)?;
  let entry = find_entry(&doc, "UNITCELL", name)?.ok_or_else(|| XmlError::UnknownEntry(name.to_string()))?;
  parse_unitcell(entry)
}

pub fn read_graph_from_str(xml: &str, name: &str) -> Result<Graph> {
  let doc = Document::parse(xml)?;
  let entry = find_entry(&doc, "GRAPH", name)?.ok_or_else(|| XmlError::UnknownEntry(name.to_string()))?;
  parse_graph(entry)
}

pub fn read_basis_from_file(path: impl AsRef<Path>, name: &str) -> Result<Basis> {
  let xml = fs::read_to_string(path)?;
  read_basis_from_str(&xml, name)
}

pub fn read_unitcell_from_file(path: impl AsRef<Path>, name: &str) -> Result<Unitcell> {
  let xml = fs::read_to_string(path)?;
  read_unitcell_from_str(&xml, name)
}

pub fn read_graph_from_file(path: impl AsRef<Path>, name: &str) -> Result<Graph> {
  let xml = fs::read_to_string(path)?;
  read_graph_from_str(&xml, name)
}

pub fn write_basis_to_string(name: &str, basis: &Basis) -> String {
  let mut xml = String::new();
  let _ = writeln!(xml, "<LATTICES>");
  let _ = writeln!(xml, "  <LATTICE name=\"{}\" dimension=\"{}\">", escape_attr(name), basis.dimension());
  let _ = writeln!(xml, "    <BASIS>");
  for column in 0..basis.dimension() {
    let mut vector = String::new();
    for row in 0..basis.dimension() {
      if row > 0 {
        vector.push(' ');
      }
      let _ = write!(vector, "{}", basis.basis_vectors()[(row, column)]);
    }
    let _ = writeln!(xml, "      <VECTOR>{}</VECTOR>", vector);
  }
  let _ = writeln!(xml, "    </BASIS>");
  let _ = writeln!(xml, "  </LATTICE>");
  let _ = writeln!(xml, "</LATTICES>");
  xml
}

pub fn write_unitcell_to_string(name: &str, cell: &Unitcell) -> String {
  let mut xml = String::new();
  let _ = writeln!(xml, "<LATTICES>");
  let _ = writeln!(xml, "  <UNITCELL name=\"{}\" dimension=\"{}\" vertices=\"{}\">", escape_attr(name), cell.dimension(), cell.num_sites());
  for site in 0..cell.num_sites() {
    let site = cell.site(site);
    let _ = write!(xml, "    <VERTEX type=\"{}\">", site.site_type);
    let _ = write!(xml, "<COORDINATE>");
    for (index, value) in site.coordinate.iter().enumerate() {
      if index > 0 {
        xml.push(' ');
      }
      let _ = write!(xml, "{}", value);
    }
    let _ = writeln!(xml, "</COORDINATE></VERTEX>");
  }
  for bond in 0..cell.num_bonds() {
    let bond = cell.bond(bond);
    let _ = write!(xml, "    <EDGE type=\"{}\">", bond.bond_type);
    let _ = write!(xml, "<SOURCE vertex=\"{}\"/>", bond.source + 1);
    let _ = write!(xml, "<TARGET vertex=\"{}\" offset=\"", bond.target + 1);
    for (index, value) in bond.target_offset.iter().enumerate() {
      if index > 0 {
        xml.push(' ');
      }
      let _ = write!(xml, "{}", value);
    }
    let _ = writeln!(xml, "\"/></EDGE>");
  }
  let _ = writeln!(xml, "  </UNITCELL>");
  let _ = writeln!(xml, "</LATTICES>");
  xml
}

pub fn write_graph_to_string(name: &str, graph: &Graph) -> String {
  let mut xml = String::new();
  let _ = writeln!(xml, "<LATTICES>");
  if graph.dimension() > 0 {
    let _ = writeln!(xml, "  <GRAPH name=\"{}\" dimension=\"{}\" vertices=\"{}\">", escape_attr(name), graph.dimension(), graph.num_sites());
  } else {
    let _ = writeln!(xml, "  <GRAPH name=\"{}\" vertices=\"{}\">", escape_attr(name), graph.num_sites());
  }
  for site in 0..graph.num_sites() {
    if graph.dimension() > 0 {
      let _ = write!(xml, "    <VERTEX type=\"{}\">", graph.site_type(site));
      let _ = write!(xml, "<COORDINATE>");
      for (index, value) in graph.coordinate(site).iter().enumerate() {
        if index > 0 {
          xml.push(' ');
        }
        let _ = write!(xml, "{}", value);
      }
      let _ = writeln!(xml, "</COORDINATE></VERTEX>");
    } else {
      let _ = writeln!(xml, "    <VERTEX type=\"{}\"/>", graph.site_type(site));
    }
  }
  for bond in 0..graph.num_bonds() {
    let _ = writeln!(xml, "    <EDGE type=\"{}\" source=\"{}\" target=\"{}\"/>", graph.bond_type(bond), graph.source(bond) + 1, graph.target(bond) + 1);
  }
  let _ = writeln!(xml, "  </GRAPH>");
  let _ = writeln!(xml, "</LATTICES>");
  xml
}

fn parse_basis(entry: Node<'_, '_>) -> Result<Basis> {
  let dimension = parse_optional_usize(entry, "dimension");
  let basis_node = entry.children().find(|child| child.is_element() && child.tag_name().name() == "BASIS")
    .ok_or(XmlError::InvalidFormat("missing BASIS element"))?;
  let mut vectors = Vec::new();
  for vector_node in basis_node.children().filter(|child| child.is_element() && child.tag_name().name() == "VECTOR") {
    vectors.push(parse_number_list::<f64>(vector_node.text().unwrap_or(""))?);
  }
  let dim = match dimension {
    Some(value) => value,
    None => vectors.len(),
  };
  if dim == 0 || vectors.len() != dim {
    return Err(XmlError::DimensionMismatch("basis dimension mismatch"));
  }
  for vector in &vectors {
    if vector.len() != dim {
      return Err(XmlError::DimensionMismatch("basis dimension mismatch"));
    }
  }
  let mut matrix = BasisMatrix::zeros(dim, dim);
  for (column, vector) in vectors.iter().enumerate() {
    for (row, value) in vector.iter().enumerate() {
      matrix[(row, column)] = *value;
    }
  }
  Ok(Basis::new(matrix))
}

fn parse_unitcell(entry: Node<'_, '_>) -> Result<Unitcell> {
  let dimension = parse_required_usize(entry, "dimension")?;
  let mut cell = Unitcell::new(dimension);
  let site_count = parse_optional_usize(entry, "vertices").unwrap_or(0);

  for vertex in entry.children().filter(|child| child.is_element() && child.tag_name().name() == "VERTEX") {
    let site_type = parse_optional_i32(vertex, "type").unwrap_or(0);
    let mut coordinate = CoordinateVector::zeros(dimension);
    let mut found_coordinate = false;
    for coordinate_node in vertex.children().filter(|child| child.is_element() && child.tag_name().name() == "COORDINATE") {
      if found_coordinate {
        return Err(XmlError::InvalidFormat("duplicated COORDINATE tag"));
      }
      let values = parse_number_list::<f64>(coordinate_node.text().unwrap_or(""))?;
      if values.len() != dimension {
        return Err(XmlError::DimensionMismatch("site coordinate dimension mismatch"));
      }
      for (index, value) in values.into_iter().enumerate() {
        coordinate[index] = value;
      }
      found_coordinate = true;
    }
    cell.add_site(coordinate, site_type);
  }

  if cell.num_sites() > 0 {
    if site_count > 0 && cell.num_sites() != site_count {
      return Err(XmlError::DimensionMismatch("inconsistent number of sites"));
    }
  } else {
    for _ in 0..site_count {
      cell.add_site(CoordinateVector::zeros(dimension), 0);
    }
  }

  for edge in entry.children().filter(|child| child.is_element() && child.tag_name().name() == "EDGE") {
    let bond_type = parse_optional_i32(edge, "type").unwrap_or(0);
    let source = parse_edge_vertex(edge, "SOURCE")?;
    let target = parse_edge_vertex(edge, "TARGET")?;
    let mut source_offset = OffsetVector::zeros(dimension);
    let mut target_offset = OffsetVector::zeros(dimension);
    for source_node in edge.children().filter(|child| child.is_element() && child.tag_name().name() == "SOURCE") {
      if let Some(vertex) = source_node.attribute("vertex") {
        let parsed = vertex.parse::<usize>().map_err(|_| XmlError::InvalidFormat("invalid SOURCE vertex"))?;
        if parsed == 0 {
          return Err(XmlError::InvalidFormat("SOURCE vertex is 1-based"));
        }
      }
      if let Some(offset) = source_node.attribute("offset") {
        let values = parse_number_list::<i64>(offset)?;
        if values.len() != dimension {
          return Err(XmlError::DimensionMismatch("SOURCE offset dimension mismatch"));
        }
        for (index, value) in values.into_iter().enumerate() {
          source_offset[index] = value;
        }
      }
    }
    for target_node in edge.children().filter(|child| child.is_element() && child.tag_name().name() == "TARGET") {
      if let Some(vertex) = target_node.attribute("vertex") {
        let parsed = vertex.parse::<usize>().map_err(|_| XmlError::InvalidFormat("invalid TARGET vertex"))?;
        if parsed == 0 {
          return Err(XmlError::InvalidFormat("TARGET vertex is 1-based"));
        }
      }
      if let Some(offset) = target_node.attribute("offset") {
        let values = parse_number_list::<i64>(offset)?;
        if values.len() != dimension {
          return Err(XmlError::DimensionMismatch("TARGET offset dimension mismatch"));
        }
        for (index, value) in values.into_iter().enumerate() {
          target_offset[index] = value;
        }
      }
    }
    let final_offset = target_offset - source_offset;
    cell.add_bond(source - 1, target - 1, final_offset, bond_type);
  }

  Ok(cell)
}

fn parse_graph(entry: Node<'_, '_>) -> Result<Graph> {
  let dimension = parse_optional_usize(entry, "dimension").unwrap_or(0);
  let mut graph = Graph::new(dimension);
  let site_count = parse_optional_usize(entry, "vertices").unwrap_or(0);

  for vertex in entry.children().filter(|child| child.is_element() && child.tag_name().name() == "VERTEX") {
    let site_type = parse_optional_i32(vertex, "type").unwrap_or(0);
    let coordinate = if dimension > 0 {
      let coordinate_node = vertex.children().find(|child| child.is_element() && child.tag_name().name() == "COORDINATE")
        .ok_or(XmlError::InvalidFormat("missing COORDINATE tag"))?;
      let values = parse_number_list::<f64>(coordinate_node.text().unwrap_or(""))?;
      if values.len() != dimension {
        return Err(XmlError::DimensionMismatch("site coordinate dimension mismatch"));
      }
      let mut coordinate = CoordinateVector::zeros(dimension);
      for (index, value) in values.into_iter().enumerate() {
        coordinate[index] = value;
      }
      coordinate
    } else {
      CoordinateVector::zeros(0)
    };
    graph.add_site(coordinate, site_type);
  }

  if graph.num_sites() > 0 {
    if site_count > 0 && graph.num_sites() != site_count {
      return Err(XmlError::DimensionMismatch("inconsistent number of sites"));
    }
  } else {
    for _ in 0..site_count {
      graph.add_site(CoordinateVector::zeros(dimension), 0);
    }
  }

  for edge in entry.children().filter(|child| child.is_element() && child.tag_name().name() == "EDGE") {
    let source = parse_required_usize_attr(edge, "source")?;
    let target = parse_required_usize_attr(edge, "target")?;
    let bond_type = parse_optional_i32(edge, "type").unwrap_or(0);
    graph.add_bond(source - 1, target - 1, bond_type);
  }

  Ok(graph)
}

fn find_entry<'a>(doc: &'a Document<'a>, tag: &str, name: &str) -> Result<Option<Node<'a, 'a>>> {
  let root = doc.root_element();
  if root.tag_name().name() != "LATTICES" {
    return Err(XmlError::InvalidFormat("root element must be LATTICES"));
  }
  Ok(root.children().find(|child| {
    child.is_element()
      && child.tag_name().name() == tag
      && child.attribute("name") == Some(name)
  }))
}

fn parse_required_usize(node: Node<'_, '_>, attr: &'static str) -> Result<usize> {
  parse_required_attr(node, attr)?.parse::<usize>().map_err(|_| XmlError::InvalidFormat(attr))
}

fn parse_required_usize_attr(node: Node<'_, '_>, attr: &'static str) -> Result<usize> {
  node.attribute(attr)
    .ok_or(XmlError::MissingAttribute(attr))?
    .parse::<usize>()
    .map_err(|_| XmlError::InvalidFormat(attr))
}

fn parse_required_attr<'a>(node: Node<'a, 'a>, attr: &'static str) -> Result<&'a str> {
  node.attribute(attr).ok_or(XmlError::MissingAttribute(attr))
}

fn parse_optional_usize(node: Node<'_, '_>, attr: &'static str) -> Option<usize> {
  node.attribute(attr).and_then(|value| value.parse::<usize>().ok())
}

fn parse_optional_i32(node: Node<'_, '_>, attr: &'static str) -> Option<i32> {
  node.attribute(attr).and_then(|value| value.parse::<i32>().ok())
}

fn parse_edge_vertex(edge: Node<'_, '_>, tag: &'static str) -> Result<usize> {
  let node = edge.children().find(|child| child.is_element() && child.tag_name().name() == tag)
    .ok_or(XmlError::MissingAttribute(tag))?;
  let vertex = node.attribute("vertex").ok_or(XmlError::MissingAttribute("vertex"))?.parse::<usize>().map_err(|_| XmlError::InvalidFormat("invalid vertex"))?;
  if vertex == 0 {
    return Err(XmlError::InvalidFormat("vertex is 1-based"));
  }
  Ok(vertex)
}

fn parse_number_list<T>(text: &str) -> Result<Vec<T>>
where
  T: std::str::FromStr,
{
  let mut values = Vec::new();
  for token in text.split_whitespace() {
    values.push(token.parse::<T>().map_err(|_| XmlError::InvalidFormat("invalid numeric value"))?);
  }
  Ok(values)
}

fn escape_attr(text: &str) -> String {
  text.replace('&', "&amp;")
    .replace('"', "&quot;")
    .replace('<', "&lt;")
    .replace('>', "&gt;")
}
