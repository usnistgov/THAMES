<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet version="1.0"
xmlns:xsl=?http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/">
  <html>
  <body>
  <h2>Microstructure Phase Chemistry and Kinetics</h2>
  <h3>General Cement Characteristics</h3>
  <table border="1">
    <tr>
      <th>Parameter</th>
      <th>Value</th>
    </tr>
    <tr>
      <td>Blaine fineness (m<sup>2</sup>/kg)</td>
      <td><xsl:value-of select:"chemistry_data/blaine"/></td>
    </tr>
    <tr>
      <td>Ref. Blaine fineness (m<sup>2</sup>/kg)</td>
      <td><xsl:value-of select:"chemistry_data/refblaine"/></td>
    </tr>
    <tr>
      <td>w/c mass ratio</td>
      <td><xsl:value-of select:"chemistry_data/wcRatio"/></td>
    </tr>
    <tr>
      <td>Temperature</td>
      <td><xsl:value-of select:"chemistry_data/temperature"/></td>
    </tr>
    <tr>
      <td>Ref. temperature</td>
      <td><xsl:value-of select:"chemistry_data/reftemperature"/></td>
    </tr>
  </table>
  <h3>Phase Characteristics</h3>
  <table border="1">
    <tr>
      <th>Phase</th>
      <th>Parameter</th>
      <th>Value</th>
    </tr>
    <xsl:for-each select="chemistry_data/phase">
    <tr>
      <td><xsl:value-of select="thamesname"/></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <td></td>
      <td>ID number</td>
      <td><xsl:value-of select:"id"/></td>
    </tr>
    <xsl:for-each select="interface_data"/>
    <tr>
      <td></td>
      <td>Random growth factor</td>
      <td><xsl:value-of select:"randomgrowth"/></td>
    </tr>
    <xsl:for-each select="growthtemplate"/>
    <tr>
      <td></td>
      <td>Growth template</td>
      <td><xsl:value-of select:"growthtemplate"/></td>
    </tr>
    </xsl:for-each>
    <xsl:for-each select="affinity"/>
    <tr>
      <td></td>
      <td>Affinity ID, Value</td>
      <td><xsl:value-of select:"affinityphaseid"/>,<xsl:value-of select:"affinityvalue"/></td>
    </tr>
    </xsl:for-each>
    </xsl:for-each>
    <tr>
      <td></td>
      <td>Porosity</td>
      <td><xsl:value-of select:"porosity"/></td>
    </tr>
    <xsl:for-each select="impurity_data"/>
    <tr>
      <td></td>
      <td>K<sub>2</sub>O Level</td>
      <td><xsl:value-of select:"k2ocoeff"/></td>
    </tr>
    <tr>
      <td></td>
      <td>Na<sub>2</sub>O Level</td>
      <td><xsl:value-of select:"na2ocoeff"/></td>
    </tr>
    <tr>
      <td></td>
      <td>MgO Level</td>
      <td><xsl:value-of select:"mgocoeff"/></td>
    </tr>
    <tr>
      <td></td>
      <td>SO<sub>3</sub> Level</td>
      <td><xsl:value-of select:"so3coeff"/></td>
    </tr>
    </xsl:for-each>
    <xsl:for-each select="kinetic_data"/>
    <tr>
      <td></td>
      <td>Kinetic type</td>
      <td><xsl:value-of select:"type"/></td>
    </tr>
    <tr>
      <td></td>
      <td>Scaled mass (g/100g)</td>
      <td><xsl:value-of select:"scaledmass"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK K<sub>1</sub> coeff</td>
      <td><xsl:value-of select:"k1"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK K<sub>2</sub> coeff</td>
      <td><xsl:value-of select:"k2"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK K<sub>3</sub> coeff</td>
      <td><xsl:value-of select:"k3"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK N<sub>1</sub> coeff</td>
      <td><xsl:value-of select:"n1"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK N<sub>3</sub> coeff</td>
      <td><xsl:value-of select:"n3"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK Activation Energy (J/mol)</td>
      <td><xsl:value-of select:"Ea"/></td>
    </tr>
    <tr>
      <td></td>
      <td>PK Critical DOH</td>
      <td><xsl:value-of select:"critdoh"/></td>
    </tr>
    <xsl:for-each select="Rd"/>
    <tr>
      <td></td>
      <td>Partition coeff for <xsl:value-of select:"Rdelement"/></td>
      <td><xsl:value-of select:"Rdvalue"/></td>
    </tr>
    </xsl:for-each>
    </xsl:for-each>
    </xsl:for-each>
  </table>
  </body>
  </html>
</xsl:template>
</xsl:stylesheet>
