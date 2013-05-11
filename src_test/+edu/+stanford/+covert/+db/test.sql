/*******************************
* test SQL procedures for MySQLDatabase_Test class
*
* Last Updated: 2/14/2012
* Author: Jonathan Karr
* Affiliation: Covert Lab, Department of Bioengineering, Stanford University
*******************************/

DROP TABLE IF EXISTS `test`;

DELIMITER $$

CREATE TABLE `test` (
  `id` bigint(20) default NULL,
  `data` longblob
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_ci$$

DROP PROCEDURE IF EXISTS `testBlobIn` $$
CREATE PROCEDURE `testBlobIn` (IN _id bigint, IN _data longblob)
BEGIN

DELETE FROM test where id=_id;

INSERT INTO test (id, data) 
VALUES (_id, _data);

END $$

DROP PROCEDURE IF EXISTS `testBlobOut` $$
CREATE PROCEDURE `testBlobOut` (IN _id bigint)
BEGIN

SELECT data
FROM test
WHERE id = _id;

END $$