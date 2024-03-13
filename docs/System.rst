######
System
######

system intro

**********
Simulation
**********

simulation intro

Information
===========

Data entry description:

* **type**: ``Simulation``, ``Information``.
* **parameters**:

  * ``name`` : ``str``, *required*.

  * ``seed`` : ``unsigned long long int``, *optional*, default: ``0xf31337Bada55D00dULL``.

* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
      "type": ["Simulation", "Information"],
      "parameters": {
        "name": "simulation name",
        "seed": 123456
      }
   }


Backup
======

Data entry description:

* **type**: ``Simulation``, ``Backup``.
* **parameters**:

  * ``backupFilePath`` : ``str``, *required*.

  * ``backupStartStep`` : ``unsigned long long int``, *optional*, default: ``0``.

  * ``backupIntervalStep`` : ``unsigned long long int``, *optional*, default: ``0``.

  * ``backupEndStep`` : ``unsigned long long int``, *optional*, default: ``MaxUnsignedLongLongInt``.

  * ``restartedFromBackup`` : ``bool``, *optional*, default: ``false``.

  * ``lastBackupStep`` : ``unsigned long long int``, *optional*.

* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
      "type": ["Simulation", "Backup"],
      "parameters": {
         "backupFilePath": "path/to/backup/file",

         "backupStartStep": 1000,
         "backupIntervalStep": 1000,
         "backupEndStep": 100000,

         "restartedFromBackup": false,
         "lastBackupStep": 0
      }
   }




********
Examples
********

.. code-block:: json

   {
      "system": {
         "information": {
            "name": "system name"
         }
      },

      //Additionnal sections
      "...":"..."
   }

.. code-block:: json

   {
      "system": {
         "information": {
            "name": "system name",
            "seed": 123456
         },
         "backup": {
            "backupFilePath": "path/to/backup/file",

            "backupStartStep": 1000,
            "backupIntervalStep": 1000,
            "backupEndStep": 100000,

            "restartedFromBackup": false,
            "lastBackupStep": 0
         }
      },

      //Additionnal sections
      "...":"..."

   }



