/**
@file valid.h
@brief Declares a function to decide if an XML file is valid compared to a schema.
*/
#ifndef VALID
#define VALID

#include <iostream>
#include <iomanip>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xmlschemas.h>

/**
@brief Decide whether an XML file is valid compared to a given schema.

@todo Make this function return a boolean instead of an int.

@param doc is a pointer to the XML file
@param schema_filename is the name of the XML schema file
@return 1 if the XML file is valid, 0 otherwise
*/
int is_xml_valid (const xmlDocPtr doc,
                  const char *schema_filename);

#endif
