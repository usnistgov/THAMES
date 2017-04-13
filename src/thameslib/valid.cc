/**
@file valid.cc
@brief Function to determine the validity of an XML file

@todo Make this return a boolean instead of an int
*/
#include "valid.h"

int is_xml_valid (const xmlDocPtr doc,
                  const char *schema_filename)
{
    xmlDocPtr schema_doc = xmlReadFile(schema_filename, NULL, XML_PARSE_NONET);
    if (schema_doc == NULL) {

        //
        // The schema cannot be loaded or is not well-formed */
        //

        std::cout << "The schema cannot be loaded or is not well-formed "
                  << std::endl;
        return -1;
    }

    xmlSchemaParserCtxtPtr parser_ctxt = xmlSchemaNewDocParserCtxt(schema_doc);
    if (parser_ctxt == NULL) {
        std::cout << "Unable to create a parser context for the schema"
                  << std::endl;
        xmlFreeDoc(schema_doc);
        return -2;
    }

    xmlSchemaPtr schema = xmlSchemaParse(parser_ctxt);
    if (schema == NULL) {

        //
        // The schema itself is not valid */
        //

        std::cout << "The schema itself is not valid" << std::endl;
        xmlSchemaFreeParserCtxt(parser_ctxt);
        xmlFreeDoc(schema_doc);
        return -3;
    }

    xmlSchemaValidCtxtPtr valid_ctxt = xmlSchemaNewValidCtxt(schema);
    if (valid_ctxt == NULL) {

        //
        // Unable to create a validation context for the schema
        //

        std::cout << "Unable to create a validation context for the schema"
                  << std::endl;
        xmlSchemaFree(schema);
        xmlSchemaFreeParserCtxt(parser_ctxt);
        xmlFreeDoc(schema_doc);
        return -4; 
    }

    int is_valid = (xmlSchemaValidateDoc(valid_ctxt, doc) == 0);
    xmlSchemaFreeValidCtxt(valid_ctxt);
    xmlSchemaFree(schema);
    xmlSchemaFreeParserCtxt(parser_ctxt);
    xmlFreeDoc(schema_doc);

    //
    // Force the return value to be non-negative on success */
    //

    std::cout << " force the return value to be non-negative on success "
              << is_valid << std::endl;

    return is_valid ? 1 : 0;
}


