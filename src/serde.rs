//! Serialization and deserialization using different formats and data types.
//!
//! # Examples
//!
//! Suppose you have a struct that you want to serialize to json and you want the BigFloat values contained in the struct to be exported as formatted strings.
//! The following is an example of how this can be achieved using `serde` field attributes and `num_bigfloat::serde::str` module:
//!
//! ``` rust
//! use num_bigfloat::BigFloat;
//! use serde::{Serialize, Deserialize};
//!
//! // A struct with a field of type BigFloat
//! #[derive(Serialize, Deserialize)]
//! struct SomeStruct {
//!
//!     #[serde(with = "num_bigfloat::serde::str")]
//!     pub f: BigFloat,
//! }
//!
//! // Value to be serialized
//! let f = BigFloat::parse("-1.234567890123456789012345678901234567890e-12").unwrap();
//! let val = SomeStruct {
//!     f,
//! };
//!
//! // Serialization
//! let json = serde_json::to_string(&val).unwrap();
//!
//! // Result
//! assert_eq!("{\"f\":\"-1.234567890123456789012345678901234567890e-12\"}", json);
//! ```

pub mod str {

    //! Serialization to string and deserialization from string.

    use crate::util::WritableBuf;
    use crate::BigFloat;
    use core::str::from_utf8_unchecked;
    use serde::Deserializer;
    use serde::Serializer;

    use serde::de::{self, Visitor};

    struct StrVisitor;

    impl<'de> Visitor<'de> for StrVisitor {
        type Value = BigFloat;

        fn expecting(&self, formatter: &mut core::fmt::Formatter) -> core::fmt::Result {
            formatter.write_str("a string representation of BigFloat")
        }

        fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            BigFloat::parse(value).ok_or_else(|| de::Error::custom("Failed to parse BigFloat"))
        }
    }

    /// Serialize `f` to a string using given serializer.
    pub fn serialize<S>(f: &BigFloat, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut buf = [0; 64];
        let mut w = WritableBuf::new(&mut buf);

        f.write_str(&mut w).unwrap();

        let written_len = w.len();

        let s = unsafe { from_utf8_unchecked(&buf[..written_len]) };

        serializer.serialize_str(s)
    }

    /// Deserialize BigFloat from a string using given deserializer.
    pub fn deserialize<'de, D>(deserializer: D) -> Result<BigFloat, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(StrVisitor)
    }
}

#[cfg(test)]
mod tests {

    use crate::BigFloat;
    use serde::{Deserialize, Serialize};

    #[derive(Serialize, Deserialize)]
    struct Stub {
        #[serde(with = "crate::serde::str")]
        pub f: BigFloat,
    }

    #[test]
    fn test_serde_str() {
        // serialize
        let f = BigFloat::parse("-1.234567890123456789012345678901234567890e-12").unwrap();
        let stub = Stub { f };

        let json = serde_json::to_string(&stub).unwrap();

        assert_eq!(
            "{\"f\":\"-1.234567890123456789012345678901234567890e-12\"}",
            json
        );

        // deserialize
        let f = BigFloat::parse("+9.123456789012345678901234567890123456789e+12").unwrap();
        let json = "{
            \"f\": \"9.123456789012345678901234567890123456789e+12\"
        }";

        let stub: Stub = serde_json::from_str(json).unwrap();

        assert!(stub.f.cmp(&f) == Some(0));
    }
}
