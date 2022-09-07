!> \file io_lib.f90
!! Input library for user keywords read from the bump file.
!///////////////////////////////////////////////////////////////////////////////////////
!////////									////////
!////////	io_lib.f90							////////
!////////									////////
!////////	contains:							////////
!////////									////////
!////////	   CHARACTER(LEN=10) FUNCTION keyword_in_(input_line, &		////////
!////////	   SUBROUTINE skip_comments(input_line,iunit,eof)		////////
!////////	   DOUBLE PRECISION FUNCTION get_dp_value(input_line)		////////
!////////	   INTEGER FUNCTION get_int_value(input_line)			////////
!////////	   CHARACTER(LEN=132) FUNCTION get_string_value(input_line)	////////
!////////	   LOGICAL FUNCTION get_logical_value(input_line)		////////
!////////	   SUBROUTINE write_file(in_unit,iout_unit)			////////
!////////									////////
!///////////////////////////////////////////////////////////////////////////////////////
!
!> @brief Search for and return the keyword "string" for a keyword
!! in the list of keywords stored in the module "keywords".
!!
!! if no keyword is found, then routine returns the error condition
!! "no keyword".
!!
!! INPUT:
!!
!!        input_line   = The input line of max length 132
!!
!!        keyword_list = the list of keywords to be consulted
!!
!!        num_keywords = the number of keywords in the list
!!
!! OUTPUT:
!!
!!        keyword_in_  = the keyword found in the input line
!!                      or "no keyword" if one is not found
!!
!!        kword        = the number of the keyword in the list
      CHARACTER(LEN=10) FUNCTION keyword_in_(input_line, &
                                          keyword_list,num_keywords,kword)
!
!.......................................................................................
!     search  for and return the keyword "string" for a keyword 
!     in the list of keywords stored in the module "keywords".
!     if no keyword is found, then routine returns the error condition
!     "no keyword"
!
!     INPUT:
!        input_line   = The input line of max length 132
!        keyword_list = the list of keywords to be consulted
!        num_keywords = the number of keywords in the list
!
!     OUTPUT:
!        keyword_in_  = the keyword found in the input line
!                      or "no keyword" if one is not found
!        kword        = the number of the keyword in the list
!
!.......................................................................................
!
      CHARACTER (LEN = 132)                         :: input_line
      INTEGER                                       :: length, ieq, loc, kword
      CHARACTER (LEN = 10), DIMENSION(num_keywords) :: keyword_list
!
      keyword_in_ = 'no keyword'
!
      ieq = INDEX(input_line,'=')
      IF(ieq == 0)     THEN
         length = LEN_TRIM(input_line)
      ELSE
         length = ieq - 1
      END IF

      DO k = 1,num_keywords
         kword = k
         len_keyword = len_trim(keyword_list(k)) - 1
         DO n = 1,length - len_keyword
            if(input_line(n:n+len_keyword) == keyword_list(k)) GO TO 100
         END DO
      END DO
      keyword_in_ = 'no keyword'
      RETURN
!
 100  continue
      keyword_in_ = keyword_list(k)

      RETURN
!      
      END FUNCTION keyword_in_
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Read in input lines, skipping over those with a "!" in the first
!! column or a blank line. Reads from unit "iunit".
!! returns the first non-comment/ non blank line in "input_line".
      SUBROUTINE skip_comments(input_line,iunit,eof)
!
!.......................................................................
!     Read in input lines, skipping over those with a "!" in the first 
!     column or a blank line. Reads from unit "iunit". 
!     returns the first non-comment/ non blank line in "input_line".
!.......................................................................
!
      CHARACTER (LEN = 132) :: input_line
      INTEGER               :: iunit
      LOGICAL               :: eof
!
      eof = .false.
      DO
         READ(iunit,'(a132)',end = 100) input_line
         IF(len(trim(input_line)) == 0)   CYCLE
         IF(input_line(1:1) /= '!')     EXIT
      END DO
      RETURN
!
 100  CONTINUE
      eof = .true.
!
      RETURN
!
      END SUBROUTINE skip_comments
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief "Read" the "value" of a double precision number declared
!! after an = sign in a input_line
   DOUBLE PRECISION FUNCTION get_dp_value(input_line)
!
!.......................................................................
!    "Read" the "value" of a double precision number declared
!     after an = sign in a input_line
!.......................................................................
!
      CHARACTER (len = 132) :: input_line
      DOUBLE PRECISION      :: value
!
      leq = INDEX(input_line,'=')
      READ(input_line(leq+1:132),*) value
      get_dp_value = value
!
      RETURN
      END FUNCTION get_dp_value
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief "Read" the "value" of an integer number declared
!! after an = sign in a input_line
      INTEGER FUNCTION get_int_value(input_line)
!
!.......................................................................
!    "Read" the "value" of an integer number declared
!     after an = sign in a input_line
!.......................................................................
!
      CHARACTER (LEN = 132) :: input_line
      INTEGER               :: value
!
      leq = INDEX(input_line,'=')
      READ(input_line(leq+1:132),*) value
      get_int_value = value
!
      RETURN
      END FUNCTION get_int_value
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Extracts the string after the "=" sign in an input file.
      CHARACTER(LEN=132) FUNCTION get_string_value(input_line)
!
!.......................................................................
!     Extracts the string after the "=" sign in an input file
!.......................................................................
!
      CHARACTER (LEN = 132) :: input_line
!
      length = LEN(TRIM(input_line))
      n = INDEX(input_line,'=')
      DO k = n+1,length
         nn = k
         IF(input_line(nn:nn).ne.' ')     EXIT
      END DO
      get_string_value = input_line(nn:length)
!
      RETURN
      END FUNCTION get_string_value
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief "Read" the "value" of an logical declared
!! after an = sign in a input_line.
      LOGICAL FUNCTION get_logical_value(input_line)
!
!.......................................................................
!    "Read" the "value" of an logical declared
!     after an = sign in a input_line
!.......................................................................
!
      CHARACTER (LEN = 132) :: input_line
      LOGICAL               :: value
!
      leq = INDEX(input_line,'=')
      READ(input_line(leq+1:132),*) value
      get_logical_value = value
!
      RETURN
      END FUNCTION get_logical_value
!
!///////////////////////////////////////////////////////////////////////
!
!> @brief Read in a file from in_unit, and write it out onto out_unit.
      SUBROUTINE write_file(in_unit,iout_unit)
!
!.......................................................................
!     read in a file from in_unit, and write it out onto out_unit
!.......................................................................
!
      CHARACTER (LEN = 132) :: input_line
!
      WRITE(6,*) '********************************************'
      WRITE(6,*) '*                                          *'
      WRITE(6,*) '*             input file:                  *'
      WRITE(6,*) '*                                          *'
      WRITE(6,*) '********************************************'
      WRITE(6,*) ' '
!
      DO
         READ(in_unit,'(a132)', END = 100) input_line
         WRITE(iout_unit,*) input_line
      END DO
!
 100  CONTINUE
!
      WRITE(6,*) ' '
!
      RETURN
      END SUBROUTINE write_file
